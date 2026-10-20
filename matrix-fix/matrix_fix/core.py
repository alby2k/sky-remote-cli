"""Correzione delle micro-imprecisioni geometriche di profili 2D (matrici di estrusione).

Il modulo non dipende da ezdxf: lavora su una lista di ``Contour`` (polilinee con
bulge, linee, archi, centri di cerchi) espressi in coordinate WCS 2D.

Idea di fondo
-------------
1. **Topologia**: i vertici più vicini di ``merge_tol`` vengono fusi in un unico
   nodo (chiude le micro-aperture fra entità, elimina i micro-segmenti).
2. **Direzioni**: ogni segmento rettilineo ha una direzione. Le direzioni quasi
   uguali (anche fra polilinee diverse) vengono raggruppate in un'unica direzione
   esatta (parallelismo), poi agganciate a 0°/90° o alla griglia angolare
   (orizzontali/verticali, 45°, ...), poi i gruppi quasi ortogonali fra loro
   vengono resi esattamente ortogonali.
3. **Tangenze**: gli archi tangenti (a meno di ``tangent_tol``) a segmenti
   rettilinei ricevono angolo al centro e direzione della corda tali da essere
   esattamente tangenti alle direzioni corrette.
4. **Allineamenti**: segmenti paralleli con offset quasi uguale diventano
   collineari; vertici con X (o Y) quasi uguale vengono allineati.
5. **Risoluzione**: con le direzioni fissate, tutti i vincoli sono *lineari*
   nelle coordinate dei nodi (es. ``n · (B - A) = 0``). Si cerca lo spostamento
   di norma minima che li soddisfa tutti (minimi quadrati vincolati), quindi la
   geometria si muove solo di quanto strettamente necessario.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import lsqr
from scipy.spatial import cKDTree


@dataclass
class Params:
    lin_tol: float = 0.02  # mm: errore lineare massimo considerato "micro-errore"
    angle_tol: float = 0.1  # gradi: errore angolare sempre corretto
    max_angle: float = 1.0  # gradi: limite assoluto (vale per segmenti corti)
    angle_step: float = 15.0  # gradi: griglia angolare (0 = solo 0°/90°)
    tangent_tol: float = 1.0  # gradi: oltre questo una giunzione è uno spigolo
    merge_tol: float | None = None  # mm: default = lin_tol
    align_vertices: bool = True
    max_move: float | None = None  # mm: default = 5 * lin_tol
    wall_step: float = 0.1  # mm: gli spessori di parete sono multipli di questo (0 = no)
    wall_tol: float = 0.049  # mm: scarto massimo di uno spessore dal valore arrotondato
    max_wall: float = 10.0  # mm: oltre questa distanza due linee non sono una parete

    @property
    def merge(self) -> float:
        return self.lin_tol if self.merge_tol is None else self.merge_tol

    @property
    def move_limit(self) -> float:
        return 5 * self.lin_tol if self.max_move is None else self.max_move


@dataclass
class Contour:
    """Una entità del disegno in coordinate WCS 2D.

    ``kind`` è "polyline", "line", "arc" (due punti, bulge sul primo) oppure
    "point" (centro di un cerchio: partecipa solo agli allineamenti).
    ``bulges[i]`` è il bulge del segmento che parte dal vertice i (convenzione DXF).
    """

    key: object
    points: list[tuple[float, float]]
    bulges: list[float]
    closed: bool = False
    kind: str = "polyline"


@dataclass
class FixedContour:
    key: object
    points: list[tuple[float, float]]
    bulges: list[float]
    closed: bool
    src_index: list[int]  # indice del vertice originale da cui proviene ogni vertice
    deleted: bool = False


@dataclass
class Report:
    merged_gaps: int = 0
    removed_micro_segments: int = 0
    straightened_bulges: int = 0
    direction_groups: int = 0
    snapped_lines: int = 0
    perpendicular_groups: int = 0
    tangent_arcs: int = 0
    collinear_groups: int = 0
    aligned_vertex_groups: int = 0
    skipped_constraints: int = 0
    walls: int = 0
    walls_changed: int = 0
    radii_kept: int = 0
    max_line_deviation_deg: float = 0.0
    max_displacement: float = 0.0
    residual: float = 0.0
    details: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def summary(self) -> str:
        rows = [
            ("Aperture chiuse (vertici fusi)", self.merged_gaps),
            ("Micro-segmenti rimossi", self.removed_micro_segments),
            ("Micro-archi raddrizzati", self.straightened_bulges),
            ("Gruppi di direzione (paralleli)", self.direction_groups),
            ("Segmenti con direzione corretta", self.snapped_lines),
            ("Gruppi resi ortogonali", self.perpendicular_groups),
            ("Archi resi tangenti", self.tangent_arcs),
            ("Gruppi di segmenti resi collineari", self.collinear_groups),
            ("Gruppi di vertici allineati", self.aligned_vertex_groups),
            ("Spessori di parete vincolati (di cui corretti)", f"{self.walls} ({self.walls_changed})"),
            ("Raggi di raccordo conservati", self.radii_kept),
            ("Vincoli facoltativi scartati", self.skipped_constraints),
            ("Errore angolare massimo trovato [°]", f"{self.max_line_deviation_deg:.6f}"),
            ("Spostamento massimo di un vertice [mm]", f"{self.max_displacement:.6f}"),
            ("Residuo dei vincoli dopo la correzione", f"{self.residual:.2e}"),
        ]
        width = max(len(r[0]) for r in rows)
        return "\n".join(f"{name:<{width}}  {value}" for name, value in rows)


class GeometryError(RuntimeError):
    pass


# ----------------------------------------------------------------------------
# utilità angolari

def _undirected_dist(a: float, b: float) -> float:
    """Distanza fra due direzioni non orientate (radianti, modulo pi)."""
    d = (a - b) % math.pi
    return min(d, math.pi - d)


def _wrap_pi(a: float) -> float:
    return (a + math.pi) % (2 * math.pi) - math.pi


def _normal(angle_deg: float) -> tuple[float, float]:
    """Normale unitaria alla direzione ``angle_deg``; esatta per i multipli di 90°."""
    k = angle_deg / 90.0
    if abs(k - round(k)) < 1e-12:
        return [(0.0, 1.0), (-1.0, 0.0)][int(round(k)) % 2]
    r = math.radians(angle_deg)
    return (-math.sin(r), math.cos(r))


class _UnionFind:
    def __init__(self, n: int):
        self.parent = list(range(n))

    def find(self, i: int) -> int:
        while self.parent[i] != i:
            self.parent[i] = self.parent[self.parent[i]]
            i = self.parent[i]
        return i

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)


def _cluster_1d(values: list[tuple[float, int]], tol: float) -> list[list[int]]:
    """Raggruppa valori scalari: un gruppo non si estende oltre ``tol`` dal primo."""
    groups, current, start = [], [], None
    for value, item in sorted(values):
        if start is None or value - start > tol:
            if len(current) > 1:
                groups.append(current)
            current, start = [], value
        current.append(item)
    if len(current) > 1:
        groups.append(current)
    return groups


# ----------------------------------------------------------------------------
# modello interno

@dataclass
class _Edge:
    a: int
    b: int
    bulge: float
    contour: int
    is_line: bool
    length: float = 0.0
    angle: float = 0.0  # direzione della corda A->B, radianti
    group: int = -1  # gruppo di direzione (solo segmenti)
    deviation: float = 0.0  # scarto dalla direzione corretta, radianti
    new_bulge: float | None = None
    chord: tuple[int, int, float, float] | None = None  # (da, a, direzione, angolo) se tangente


@dataclass
class _DirGroup:
    angle: float  # gradi in [0, 180)
    sx: float = 0.0
    sy: float = 0.0
    members: list[int] = field(default_factory=list)
    lmax: float = 0.0
    total: float = 0.0
    fixed: bool = False


class _Fixer:
    def __init__(self, contours: list[Contour], params: Params):
        self.contours = contours
        self.p = params
        self.report = Report()

    # -- tolleranza combinata angolo / spostamento laterale ----------------
    def _within(self, d: float, length: float) -> bool:
        if d <= math.radians(self.p.angle_tol):
            return True
        return d <= math.radians(self.p.max_angle) and length * math.sin(d) <= self.p.lin_tol

    # -- 1. topologia -------------------------------------------------------
    def build_nodes(self) -> None:
        raw, owner = [], []
        for ci, c in enumerate(self.contours):
            for vi, pt in enumerate(c.points):
                raw.append(pt)
                owner.append((ci, vi))
        self.raw = np.array(raw, dtype=float).reshape(-1, 2)
        uf = _UnionFind(len(raw))
        mergeable = [self.contours[ci].kind != "point" for ci, _ in owner]
        if len(raw) > 1:
            for i, j in cKDTree(self.raw).query_pairs(self.p.merge):
                if mergeable[i] and mergeable[j]:
                    uf.union(i, j)
        roots: dict[int, int] = {}
        self.raw_node = []
        for i in range(len(raw)):
            self.raw_node.append(roots.setdefault(uf.find(i), len(roots)))
        n = len(roots)
        sums = np.zeros((n, 2))
        counts = np.zeros(n)
        for i, node in enumerate(self.raw_node):
            sums[node] += self.raw[i]
            counts[node] += 1
        self.nodes = sums / counts[:, None]
        self.raw_index = {}
        k = 0
        for ci, c in enumerate(self.contours):
            for vi in range(len(c.points)):
                self.raw_index[(ci, vi)] = k
                k += 1

    def build_edges(self) -> None:
        self.edges: list[_Edge] = []
        self.sequences = []  # per contorno: (nodi, bulge, src_index, closed, edge ids)
        for ci, c in enumerate(self.contours):
            nodes = [self.raw_node[self.raw_index[(ci, vi)]] for vi in range(len(c.points))]
            bulges = list(c.bulges) + [0.0] * (len(nodes) - len(c.bulges))
            src = list(range(len(nodes)))
            closed = c.closed
            if c.kind == "point":
                self.sequences.append((nodes, bulges, src, False, []))
                continue
            changed = True
            while changed and len(nodes) > 1:
                changed = False
                for i in range(len(nodes) - 1):
                    if nodes[i] == nodes[i + 1]:
                        del nodes[i], bulges[i], src[i]
                        self.report.removed_micro_segments += 1
                        changed = True
                        break
                if c.kind == "polyline" and len(nodes) > 1 and nodes[0] == nodes[-1]:
                    if closed:
                        self.report.removed_micro_segments += 1
                    else:
                        closed = True  # polilinea "chiusa a vista" ma aperta: la chiudo
                    del nodes[-1], bulges[-1], src[-1]
                    changed = True
            if len(nodes) < 2:
                self.sequences.append((nodes, bulges, src, closed, None))
                continue
            edge_ids = []
            nseg = len(nodes) if closed else len(nodes) - 1
            for i in range(nseg):
                a, b = nodes[i], nodes[(i + 1) % len(nodes)]
                bulge = bulges[i]
                va = self.nodes[a]
                vb = self.nodes[b]
                chord = float(np.hypot(*(vb - va)))
                if bulge and c.kind != "arc" and abs(bulge) * chord / 2 < self.p.lin_tol / 2:
                    bulge = 0.0
                    bulges[i] = 0.0
                    self.report.straightened_bulges += 1
                e = _Edge(a, b, bulge, ci, is_line=(bulge == 0.0), length=chord,
                          angle=math.atan2(vb[1] - va[1], vb[0] - va[0]))
                edge_ids.append(len(self.edges))
                self.edges.append(e)
            self.sequences.append((nodes, bulges, src, closed, edge_ids))
        # aperture chiuse: nodi nati da punti non coincidenti
        gaps = {self.raw_node[k] for k in range(len(self.raw_node))
                if np.any(self.raw[k] != self.nodes[self.raw_node[k]])}
        self.report.merged_gaps = len(gaps)

    # -- 2. direzioni -------------------------------------------------------
    def cluster_directions(self) -> None:
        lines = sorted((i for i, e in enumerate(self.edges) if e.is_line),
                       key=lambda i: -self.edges[i].length)
        groups: list[_DirGroup] = []
        for i in lines:
            e = self.edges[i]
            theta = e.angle % math.pi
            best, best_d = None, None
            for gi, g in enumerate(groups):
                d = _undirected_dist(theta, math.radians(g.angle))
                if best_d is None or d < best_d:
                    best, best_d = gi, d
            if best is None or not self._within(best_d, e.length):
                groups.append(_DirGroup(angle=math.degrees(theta) % 180.0))
                best = len(groups) - 1
            g = groups[best]
            g.members.append(i)
            g.sx += e.length * math.cos(2 * theta)
            g.sy += e.length * math.sin(2 * theta)
            g.angle = math.degrees(math.atan2(g.sy, g.sx) / 2) % 180.0
            g.lmax = max(g.lmax, e.length)
            g.total += e.length
            e.group = best

        # angolo di progetto: se un segmento del gruppo ha già un angolo "tondo"
        # esatto (multiplo di 0.01°), quello è il valore voluto, non la media
        for g in groups:
            for i in g.members:  # già ordinati dal più lungo
                deg = math.degrees(self.edges[i].angle) % 180.0
                exact = round(deg * 100) / 100
                if abs(deg - exact) < 1e-7:
                    g.angle = exact % 180.0
                    g.fixed = True
                    break

        # aggancio alla griglia angolare (0° e 90° sempre)
        step = self.p.angle_step
        for g in groups:
            candidates = [0.0, 90.0, 180.0]
            if step > 0:
                k = round(g.angle / step)
                candidates += [(k - 1) * step, k * step, (k + 1) * step]
            target = min(candidates, key=lambda c: abs(g.angle - c))
            d = math.radians(abs(g.angle - target))
            if self._within(d, g.lmax):
                g.angle = target % 180.0
                g.fixed = True

        # ortogonalità fra gruppi non agganciati alla griglia
        for g in sorted(groups, key=lambda g: -g.total):
            if g.fixed:
                continue
            for ref in sorted((r for r in groups if r.fixed), key=lambda r: -r.total):
                target = (ref.angle + 90.0) % 180.0
                d = _undirected_dist(math.radians(g.angle), math.radians(target))
                if self._within(d, g.lmax):
                    g.angle = target
                    self.report.perpendicular_groups += 1
                    break
            g.fixed = True

        self.groups = groups
        self.report.direction_groups = len({round(g.angle, 9) for g in groups})
        for e in (self.edges[i] for i in lines):
            d = _undirected_dist(e.angle, math.radians(groups[e.group].angle))
            e.deviation = d
            if d > 1e-12:
                self.report.snapped_lines += 1
                self.report.max_line_deviation_deg = max(
                    self.report.max_line_deviation_deg, math.degrees(d))
                if d > math.radians(1e-6):
                    self.report.details.append(
                        f"entità {self.contours[e.contour].key}: segmento lungo {e.length:.4f} mm "
                        f"a {math.degrees(e.angle) % 180:.6f}° -> {groups[e.group].angle:.6f}° "
                        f"(errore {math.degrees(d):.6f}°)")

    # -- 3. tangenze --------------------------------------------------------
    def _alpha(self, e: _Edge) -> float:
        return 4 * math.atan(e.bulge)

    def _out_tangent(self, ei: int, end: int) -> float:
        """Direzione uscente dal nodo lungo l'entità (radianti)."""
        e = self.edges[ei]
        half = 0.0 if e.is_line else self._alpha(e) / 2
        if end == 0:
            return e.angle - half
        return e.angle + half + math.pi

    def tangencies(self) -> None:
        ends: dict[int, list[tuple[int, int]]] = {}
        for ei, e in enumerate(self.edges):
            ends.setdefault(e.a, []).append((ei, 0))
            ends.setdefault(e.b, []).append((ei, 1))
        link: dict[tuple[int, int], tuple[int, int]] = {}
        tol = math.radians(self.p.tangent_tol)
        for node, lst in ends.items():
            if len(lst) != 2 or lst[0][0] == lst[1][0]:
                continue
            (e1, s1), (e2, s2) = lst
            if self.edges[e1].is_line and self.edges[e2].is_line:
                continue
            diff = abs(_wrap_pi(self._out_tangent(e1, s1) - self._out_tangent(e2, s2) - math.pi))
            if diff <= tol:
                link[(e1, s1)] = (e2, s2)
                link[(e2, s2)] = (e1, s1)

        visited: set[int] = set()
        for start in range(len(self.edges)):
            if start in visited or ((start, 0) not in link and (start, 1) not in link):
                continue
            # risale all'inizio della catena (o riconosce un anello)
            cur, rev, cyclic = start, False, False
            while True:
                entry = 1 if rev else 0
                prev = link.get((cur, entry))
                if prev is None:
                    break
                cur, rev = prev[0], prev[1] == 0
                if cur == start:
                    cyclic = True
                    break
            chain = []
            while True:
                chain.append((cur, rev))
                visited.add(cur)
                nxt = link.get((cur, 0 if rev else 1))
                if nxt is None:
                    break
                cur, rev = nxt[0], nxt[1] == 1
                if cur == chain[0][0]:
                    cyclic = True
                    break
            self._solve_chain(chain, cyclic)

    def _travel_dir(self, ei: int, rev: bool) -> float:
        e = self.edges[ei]
        travel = e.angle + (math.pi if rev else 0.0)
        a = math.radians(self.groups[e.group].angle)
        return a if abs(_wrap_pi(travel - a)) <= math.pi / 2 else a + math.pi

    def _solve_chain(self, chain: list[tuple[int, bool]], cyclic: bool) -> None:
        anchors = [k for k, (ei, _) in enumerate(chain) if self.edges[ei].is_line]
        if not anchors:
            return
        alpha0 = {k: self._alpha(self.edges[ei]) * (-1 if rev else 1)
                  for k, (ei, rev) in enumerate(chain) if not self.edges[ei].is_line}
        result: dict[int, tuple[float, float]] = {}  # k -> (alpha, psi) nel verso di percorrenza
        n = len(chain)

        spans = list(zip(anchors, anchors[1:]))
        if cyclic:
            spans.append((anchors[-1], anchors[0] + n))
        for i, j in spans:
            arcs = [k % n for k in range(i + 1, j)]
            if not arcs:
                continue
            t_i = self._travel_dir(*chain[i])
            t_j = self._travel_dir(*chain[j % n])
            s0 = sum(alpha0[k] for k in arcs)
            delta = s0 + _wrap_pi(t_j - t_i - s0)
            corr = delta - s0
            if abs(corr) > 2 * math.radians(self.p.tangent_tol):
                continue
            weight = sum(abs(alpha0[k]) for k in arcs) or 1.0
            t = t_i
            for k in arcs:
                alpha = alpha0[k] + corr * abs(alpha0[k]) / weight
                result[k] = (alpha, t + alpha / 2)
                t += alpha
        if not cyclic:
            t = self._travel_dir(*chain[anchors[-1]])
            for k in range(anchors[-1] + 1, n):
                result[k] = (alpha0[k], t + alpha0[k] / 2)
                t += alpha0[k]
            t = self._travel_dir(*chain[anchors[0]])
            for k in range(anchors[0] - 1, -1, -1):
                result[k] = (alpha0[k], t - alpha0[k] / 2)
                t -= alpha0[k]

        for k, (alpha, psi) in result.items():
            ei, rev = chain[k]
            e = self.edges[ei]
            e.new_bulge = math.tan(alpha / 4) * (-1 if rev else 1)
            e.chord = (e.b, e.a, psi, alpha) if rev else (e.a, e.b, psi, alpha)
            self.report.tangent_arcs += 1

    # -- 4/5. vincoli e risoluzione -----------------------------------------
    def _material_test(self):
        """Funzione punto -> True se il punto è nel pieno del profilo (regola pari/dispari
        sui contorni chiusi). Senza contorni chiusi non si può sapere: sempre True."""
        polys = []
        for nodes, bulges, _src, closed, edge_ids in self.sequences:
            if not closed or not edge_ids:
                continue
            pts = []
            for k, ei in enumerate(edge_ids):
                e = self.edges[ei]
                pa = self.nodes[e.a]
                pts.append(pa)
                if not e.is_line:
                    pb = self.nodes[e.b]
                    chord = float(np.hypot(*(pb - pa)))
                    u = (pb - pa) / chord
                    b = e.bulge
                    h = chord * (1 - b * b) / (4 * b)
                    c = (pa + pb) / 2 + np.array([-u[1], u[0]]) * h
                    a0 = math.atan2(pa[1] - c[1], pa[0] - c[0])
                    alpha = 4 * math.atan(b)
                    rad = float(np.hypot(*(pa - c)))
                    for t in np.linspace(0, 1, 17)[1:-1]:
                        a = a0 + alpha * t
                        pts.append(c + rad * np.array([math.cos(a), math.sin(a)]))
            polys.append(np.array(pts))
        if not polys:
            return lambda pt: True

        def inside(pt) -> bool:
            count = 0
            for poly in polys:
                x, y = poly[:, 0], poly[:, 1]
                x2, y2 = np.roll(x, -1), np.roll(y, -1)
                cross = (y > pt[1]) != (y2 > pt[1])
                with np.errstate(divide="ignore", invalid="ignore"):
                    xi = x + (pt[1] - y) * (x2 - x) / (y2 - y)
                count += int(np.count_nonzero(cross & (pt[0] < xi)) % 2)
            return count % 2 == 1
        return inside

    def build_constraints(self) -> None:
        # ogni vincolo è (coefficienti, termine noto) con un'etichetta: None = vincolo
        # essenziale (direzioni, tangenze), altrimenti descrive un vincolo facoltativo,
        # che viene scartato se incompatibile; l'ordine delle liste è la priorità
        self.radius_tags: dict[str, tuple[_Edge, float]] = {}
        core: list = []
        walls: list = []
        radii: list = []
        collinear: list = []
        aligned: list = []

        def diff_row(target, normal, a, b, rhs=0.0, tag=None):
            nx, ny = normal
            row = []
            if nx:
                row += [(2 * b, nx), (2 * a, -nx)]
            if ny:
                row += [(2 * b + 1, ny), (2 * a + 1, -ny)]
            target.append((row, rhs, tag))

        def where(node):
            return f"({self.nodes[node][0]:.4f}, {self.nodes[node][1]:.4f})"

        for e in self.edges:
            if e.is_line:
                diff_row(core, _normal(self.groups[e.group].angle), e.a, e.b)
            elif e.chord is not None:
                p, q, psi, _alpha = e.chord
                diff_row(core, (-math.sin(psi), math.cos(psi)), p, q)

        # raggi: un arco con corda di direzione nota ha corda lineare nei nodi, quindi
        # il raggio di progetto (arrotondato al centesimo se ci è vicino) si conserva
        for e in self.edges:
            if e.chord is None:
                continue
            p, q, psi, alpha = e.chord
            r0 = e.length / (2 * abs(math.sin(2 * math.atan(e.bulge))))
            r = round(r0, 2) if abs(r0 - round(r0, 2)) <= 1e-3 else r0
            if abs(r - r0) > 1e-9:
                self.report.details.append(f"raggio {r0:.6f} -> {r:.2f} vicino a {where(p)}")
            self.report.radii_kept += 1
            self.radius_tags[f"raggio R{r:g} vicino a {where(p)}"] = (e, r)
            radii.append(([(2 * q, math.cos(psi)), (2 * p, -math.cos(psi)),
                           (2 * q + 1, math.sin(psi)), (2 * p + 1, -math.sin(psi))],
                          2 * r * abs(math.sin(alpha / 2)), f"raggio R{r:g} vicino a {where(p)}"))

        by_angle: dict[float, list[int]] = {}
        for i, e in enumerate(self.edges):
            if e.is_line:
                by_angle.setdefault(round(self.groups[e.group].angle, 9), []).append(i)
        in_material = self._material_test() if self.p.wall_step > 0 else None
        for angle, members in by_angle.items():
            n = _normal(angle)
            u = (n[1], -n[0])
            offset, span = {}, {}
            for i in members:
                e = self.edges[i]
                pa, pb = self.nodes[e.a], self.nodes[e.b]
                offset[i] = n[0] * (pa[0] + pb[0]) / 2 + n[1] * (pa[1] + pb[1]) / 2
                sa, sb = u[0] * pa[0] + u[1] * pa[1], u[0] * pb[0] + u[1] * pb[1]
                span[i] = (min(sa, sb), max(sa, sb))

            # collinearità: segmenti con la stessa direzione e offset quasi uguale
            uf = _UnionFind(len(self.edges))
            for grp in _cluster_1d([(offset[i], i) for i in members], self.p.lin_tol):
                ref = self.edges[grp[0]].a
                nodes_in = {self.edges[i].a for i in grp} | {self.edges[i].b for i in grp}
                if len(nodes_in) > 2:
                    self.report.collinear_groups += 1
                tag = f"collinearità a {angle:g}° vicino a {where(ref)}"
                for i in grp[1:]:
                    uf.union(grp[0], i)
                    if self.edges[i].a != ref:
                        diff_row(collinear, n, ref, self.edges[i].a, tag=tag)

            # spessori di parete: due segmenti paralleli affacciati, con il pieno in
            # mezzo e nessun altro segmento fra loro, distano un multiplo di wall_step
            if in_material is None:
                continue
            step = self.p.wall_step
            ordered = sorted(members, key=lambda i: offset[i])
            pairs = []
            for x, i in enumerate(ordered):
                for y in range(x + 1, len(ordered)):
                    j = ordered[y]
                    d = offset[j] - offset[i]
                    if d > self.p.max_wall:
                        break
                    if d <= self.p.lin_tol:
                        continue
                    lo, hi = max(span[i][0], span[j][0]), min(span[i][1], span[j][1])
                    if hi - lo <= self.p.lin_tol:
                        continue
                    t = round(round(d / step) * step, 9)
                    if t <= 0 or abs(d - t) > self.p.wall_tol:
                        continue
                    if any(offset[i] + self.p.lin_tol < offset[k] < offset[j] - self.p.lin_tol
                           and min(hi, span[k][1]) - max(lo, span[k][0]) > 0
                           for k in ordered[x + 1:y]):
                        continue  # c'è un altro segmento in mezzo: non è una parete
                    s_mid, o_mid = (lo + hi) / 2, (offset[i] + offset[j]) / 2
                    if not in_material((u[0] * s_mid + n[0] * o_mid, u[1] * s_mid + n[1] * o_mid)):
                        continue  # in mezzo c'è il vuoto (cava, luce): non è una parete
                    pairs.append((d, i, j, t))
            for d, i, j, t in sorted(pairs):
                if uf.find(i) == uf.find(j):
                    continue  # già legati: un secondo vincolo potrebbe contraddire il primo
                uf.union(i, j)
                a = self.edges[i].a
                self.report.walls += 1
                if abs(d - t) > 1e-9:
                    self.report.walls_changed += 1
                    self.report.details.append(
                        f"spessore parete a {angle:g}° vicino a {where(a)}: {d:.6f} -> {t:g} mm")
                diff_row(walls, n, a, self.edges[j].a, rhs=t,
                         tag=f"spessore {t:g} mm a {angle:g}° vicino a {where(a)}")

        # allineamento verticale/orizzontale di spigoli fra segmenti e centri di fori;
        # gli estremi degli archi sono esclusi: la loro posizione la decide la tangenza
        if self.p.align_vertices:
            adjacent = {(min(e.a, e.b), max(e.a, e.b)) for e in self.edges}
            on_arc = {n for e in self.edges if not e.is_line for n in (e.a, e.b)}
            used = sorted(({e.a for e in self.edges} | {e.b for e in self.edges}
                           | {s[0][0] for s in self.sequences if s[4] == [] and s[0]}) - on_arc)
            for axis in (0, 1):
                values = [(self.nodes[i][axis], i) for i in used]
                for grp in _cluster_1d(values, self.p.lin_tol):
                    spread = max(self.nodes[i][axis] for i in grp) - min(self.nodes[i][axis] for i in grp)
                    if spread > 0:
                        self.report.aligned_vertex_groups += 1
                    tag = f"allineamento in {'XY'[axis]} di {', '.join(where(i) for i in grp)}"
                    linked = [grp[0]]
                    for i in grp[1:]:
                        ref = next((r for r in linked if (min(r, i), max(r, i)) not in adjacent), None)
                        if ref is not None:
                            aligned.append(([(2 * i + axis, 1.0), (2 * ref + axis, -1.0)], 0.0, tag))
                            linked.append(i)
        self.constraints = core + walls + radii + collinear + aligned

    def _unchanged(self, e: _Edge) -> bool:
        if e.is_line:
            return e.deviation < 1e-11
        if e.chord is None:
            return True
        p, q, psi, _alpha = e.chord
        orig = math.atan2(*(self.nodes[q] - self.nodes[p])[::-1])
        return abs(e.new_bulge - e.bulge) < 1e-12 and abs(_wrap_pi(psi - orig)) < 1e-11

    def _max_move(self, v: np.ndarray) -> float:
        if not len(self.raw):
            return 0.0
        return float(np.hypot(*(v.reshape(-1, 2)[self.raw_node] - self.raw).T).max())

    def solve(self, cons: list) -> tuple[np.ndarray, float]:
        """Spostamento minimo che soddisfa i vincoli; restituisce (coordinate, residuo)."""
        nvar = 2 * len(self.nodes)
        v = self.nodes.reshape(-1).copy()
        if not cons:
            return v, 0.0
        data, ri, ci = [], [], []
        for r, (row, _rhs, _tag) in enumerate(cons):
            for c, val in row:
                ri.append(r)
                ci.append(c)
                data.append(val)
        b = np.array([rhs for _row, rhs, _tag in cons])
        # i vertici in cui arrivano solo entità rimaste invariate (segmenti già nella
        # direzione esatta, archi con angolo e corda invariati) sono "fidati":
        # spostarli costa molto di più, così la correzione si scarica sulle parti
        # davvero sbagliate (minimi quadrati pesati: colonne scalate di 1/sqrt(peso))
        scale = np.ones(nvar)
        exact: dict[int, bool] = {}
        for e in self.edges:
            for n in (e.a, e.b):
                exact[n] = exact.get(n, True) and self._unchanged(e)
        for n, ok in exact.items():
            if ok:
                scale[2 * n:2 * n + 2] = 1e-4
        A = csr_matrix((data, (ri, ci)), shape=(len(cons), nvar))
        As = A @ diags(scale)
        for _ in range(4):
            r = b - A @ v
            if np.max(np.abs(r)) < 1e-13:
                break
            dv = scale * lsqr(As, r, atol=1e-16, btol=1e-16, conlim=1e12,
                              iter_lim=50 * nvar + 1000)[0]
            v = v + dv
        # ciò che si è mosso meno di 1e-9 mm torna identico all'originale (bit per bit)
        v0 = self.nodes.reshape(-1)
        still = np.abs(v - v0) < 1e-9
        v[still] = v0[still]

        # uguaglianze esatte per gli assi (orizzontali/verticali "bit per bit")
        uf = _UnionFind(nvar)
        for row, rhs, _tag in cons:
            if (rhs == 0.0 and len(row) == 2 and row[0][1] == -row[1][1]
                    and row[0][0] % 2 == row[1][0] % 2):
                uf.union(row[0][0], row[1][0])
        groups: dict[int, list[int]] = {}
        for i in range(nvar):
            groups.setdefault(uf.find(i), []).append(i)
        for members in groups.values():
            if len(members) > 1:
                v[members] = float(np.mean(v[members]))
        return v, float(np.max(np.abs(b - A @ v)))

    def solve_with_fallback(self) -> np.ndarray:
        """Risolve tutto; un vincolo facoltativo che costringe a spostamenti oltre
        ``max_move`` o che contraddice gli altri viene scartato (ricerca per
        bisezione, in ordine di priorità) e segnalato."""
        limit = self.p.move_limit

        def ok(v, res):
            return self._max_move(v) <= limit and res <= 1e-9

        v, res = self.solve(self.constraints)
        if ok(v, res):
            self.report.residual = res
            return v
        core = [c for c in self.constraints if c[2] is None]
        optional: dict[str, list] = {}
        for c in self.constraints:
            if c[2] is not None:
                optional.setdefault(c[2], []).append(c)
        v, res = self.solve(core)
        if not ok(v, res) or not optional:
            self.report.residual = res
            return v
        kept: list = []
        dropped: list = []

        def try_add(batch: list[str]) -> None:
            nonlocal v, res
            if not batch:
                return
            cand = kept + [c for t in batch for c in optional[t]]
            v2, res2 = self.solve(core + cand)
            if ok(v2, res2):
                kept[:] = cand
                v, res = v2, res2
            elif len(batch) == 1:
                dropped.append((batch[0], res2, self._max_move(v2)))
            else:
                try_add(batch[: len(batch) // 2])
                try_add(batch[len(batch) // 2:])

        try_add(list(optional))
        nodes = v.reshape(-1, 2)
        for tag, res2, move in dropped:
            if tag in self.radius_tags:
                # un raggio lasciato libero per far chiudere una catena di raccordi
                # tutti a raggio fisso: se resta quasi uguale è solo arrotondamento
                e, r = self.radius_tags[tag]
                bulge = e.bulge if e.new_bulge is None else e.new_bulge
                chord = float(np.hypot(*(nodes[e.b] - nodes[e.a])))
                r_fin = chord / (2 * abs(math.sin(2 * math.atan(bulge))))
                if abs(r_fin - r) < 1e-4:
                    self.report.details.append(
                        f"{tag}: lasciato a {r_fin:.7f} per chiudere la catena di raccordi "
                        f"(scarto {abs(r_fin - r):.1e} mm)")
                    continue
            why = (f"richiederebbe spostamenti di {move:.4f} mm" if res2 <= 1e-9
                   else "in conflitto con altri vincoli")
            self.report.notes.append(f"vincolo facoltativo non applicato ({tag}): {why}")
        self.report.skipped_constraints = len(dropped)
        self.report.residual = res
        return v

    # -- esecuzione completa ------------------------------------------------
    def run(self) -> list[FixedContour]:
        self.build_nodes()
        self.build_edges()
        self.cluster_directions()
        self.tangencies()
        self.build_constraints()
        v = self.solve_with_fallback().reshape(-1, 2)

        moves = np.hypot(*(v[self.raw_node] - self.raw).T) if len(self.raw) else np.zeros(0)
        self.report.max_displacement = float(moves.max()) if len(moves) else 0.0
        if self.report.residual > 1e-7:
            self.report.warnings.append(
                f"vincoli non completamente soddisfatti (residuo {self.report.residual:.2e})")
        if self.report.max_displacement > self.p.move_limit:
            worst = int(np.argmax(moves))
            self.report.warnings.append(
                f"un vertice si sposterebbe di {self.report.max_displacement:.4f} mm "
                f"(limite {self.p.move_limit} mm), vicino a ({self.raw[worst][0]:.4f}, "
                f"{self.raw[worst][1]:.4f}): probabile vincolo errato, controllare il disegno")

        out = []
        for ci, (nodes, bulges, src, closed, edge_ids) in enumerate(self.sequences):
            c = self.contours[ci]
            if edge_ids is None:
                out.append(FixedContour(c.key, [], [], closed, [], deleted=True))
                continue
            new_bulges = list(bulges)
            for k, ei in enumerate(edge_ids):
                e = self.edges[ei]
                new_bulges[k] = 0.0 if e.is_line else (
                    e.new_bulge if e.new_bulge is not None else e.bulge)
            pts = [(float(v[n][0]), float(v[n][1])) for n in nodes]
            out.append(FixedContour(c.key, pts, new_bulges, closed, src))
        return out


def fix_contours(contours: list[Contour], params: Params | None = None
                 ) -> tuple[list[FixedContour], Report]:
    """Corregge i contorni e restituisce (contorni corretti, report).

    I contorni restituiti sono nello stesso ordine di quelli in ingresso.
    """
    fixer = _Fixer(contours, params or Params())
    fixed = fixer.run()
    return fixed, fixer.report
