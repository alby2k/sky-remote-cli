import math
import random

import ezdxf
import pytest

from matrix_fix import Contour, Params, fix_contours
from matrix_fix import dxf_io
from matrix_fix.cli import main


def fillet_polygon(corners, radius):
    """Vertici e bulge di un poligono chiuso con raccordi tangenti di raggio ``radius``."""
    pts, bulges = [], []
    n = len(corners)
    for i in range(n):
        (px, py), (cx, cy), (nx, ny) = corners[i - 1], corners[i], corners[(i + 1) % n]
        a1 = math.atan2(cy - py, cx - px)
        a2 = math.atan2(ny - cy, nx - cx)
        turn = (a2 - a1 + math.pi) % (2 * math.pi) - math.pi
        t = radius * math.tan(abs(turn) / 2)
        pts.append((cx - t * math.cos(a1), cy - t * math.sin(a1)))
        bulges.append(math.tan(turn / 4))
        pts.append((cx + t * math.cos(a2), cy + t * math.sin(a2)))
        bulges.append(0.0)
    return pts, bulges


def jitter(points, amount, rng):
    return [(x + rng.uniform(-amount, amount), y + rng.uniform(-amount, amount)) for x, y in points]


def die_profile():
    """Profilo tipo facciata: esterno con smusso a 45°, cava interna, parete a 30°."""
    outer = [(0, 0), (60, 0), (60, 30), (50, 40), (0, 40)]
    cave = [(2, 2), (58, 2), (58, 29.171573), (49.171573, 38), (2, 38)]
    s30, c30 = math.sin(math.radians(30)), math.cos(math.radians(30))
    wall = [(70, 0), (70 + 20 * c30, 20 * s30),
            (70 + 20 * c30 - 2 * s30, 20 * s30 + 2 * c30), (70 - 2 * s30, 2 * c30)]
    return outer, cave, wall


def segment_angles(points, bulges, closed):
    n = len(points)
    for i in range(n if closed else n - 1):
        if bulges[i] == 0:
            (x0, y0), (x1, y1) = points[i], points[(i + 1) % n]
            yield x1 - x0, y1 - y0


def tangent_errors(points, bulges):
    """Scarto angolare fra le tangenti in ogni vertice di una polilinea chiusa."""
    n = len(points)
    out_dir, in_dir = [], []
    for i in range(n):
        (x0, y0), (x1, y1) = points[i], points[(i + 1) % n]
        chord = math.atan2(y1 - y0, x1 - x0)
        half = 2 * math.atan(bulges[i])
        out_dir.append(chord - half)
        in_dir.append(chord + half)
    for i in range(n):
        d = (out_dir[i] - in_dir[i - 1] + math.pi) % (2 * math.pi) - math.pi
        yield abs(d), bulges[i] != 0 or bulges[i - 1] != 0


@pytest.fixture
def dirty_contours():
    rng = random.Random(42)
    outer, cave, wall = die_profile()
    o_pts, o_b = fillet_polygon(outer, 2.0)
    c_pts, c_b = fillet_polygon(cave, 0.5)
    return [
        Contour("outer", jitter(o_pts, 0.005, rng), o_b, closed=True),
        Contour("cave", jitter(c_pts, 0.005, rng), c_b, closed=True),
        Contour("wall", jitter(wall, 0.004, rng), [0.0] * 4, closed=True),
    ]


def test_directions_are_exact(dirty_contours):
    fixed, report = fix_contours(dirty_contours)
    assert not report.warnings
    assert report.max_displacement < 0.02
    allowed = [0, 30, 45, 90, 120, 135]
    for fc in fixed:
        for dx, dy in segment_angles(fc.points, fc.bulges, fc.closed):
            angle = math.degrees(math.atan2(dy, dx)) % 180
            err = min(abs(angle - a) for a in allowed + [180])
            assert err < 1e-9
            # orizzontali/verticali esatte bit per bit
            if abs(angle - 90) < 1e-6:
                assert dx == 0.0
            elif angle < 1e-6 or angle > 180 - 1e-6:
                assert dy == 0.0


def test_fillets_are_tangent(dirty_contours):
    before = max(e for c in dirty_contours[:2]
                 for e, arc in tangent_errors(c.points, c.bulges) if arc)
    assert before > 1e-5  # il jitter ha rotto le tangenze
    fixed, _ = fix_contours(dirty_contours)
    for fc in fixed[:2]:
        for err, arc in tangent_errors(fc.points, fc.bulges):
            if arc:
                assert err < 1e-9


def test_wall_thickness_lines_are_parallel(dirty_contours):
    fixed, _ = fix_contours(dirty_contours)
    wall = fixed[2]
    dirs = [math.atan2(dy, dx) % math.pi for dx, dy in segment_angles(wall.points, wall.bulges, True)]
    assert abs(dirs[0] - dirs[2]) < 1e-12
    assert abs(dirs[1] - dirs[3]) < 1e-12
    assert abs(abs(dirs[0] - dirs[1]) - math.pi / 2) < 1e-12


def test_gap_closing_and_micro_segments():
    contours = [
        # contorno fatto di LINE e ARC separati con aperture di 8 micron
        Contour("l1", [(0, 0), (10, 0.003)], [0, 0], kind="line"),
        Contour("l2", [(10.008, 0), (10, 10)], [0, 0], kind="line"),
        Contour("l3", [(10, 10), (0.004, 10.002)], [0, 0], kind="line"),
        Contour("l4", [(0, 10), (0, 0.006)], [0, 0], kind="line"),
        # polilinea aperta con fine coincidente all'inizio e micro-segmento
        Contour("pl", [(20, 0), (30, 0), (30, 0.001), (30, 10), (20, 10), (20.003, 0)],
                [0] * 6, closed=False),
    ]
    fixed, report = fix_contours(contours)
    assert not report.warnings
    lines = fixed[:4]
    for a, b in zip(lines, lines[1:] + lines[:1]):
        assert a.points[1] == b.points[0]
    pl = fixed[4]
    assert pl.closed and len(pl.points) == 4
    assert report.removed_micro_segments >= 1


def test_vertex_alignment_and_collinearity():
    contours = [
        Contour("a", [(0, 0), (10, 0), (10, 5), (0, 5)], [0] * 4, closed=True),
        Contour("b", [(20, 0.007), (30, 0.007), (30, 5), (20, 5)], [0] * 4, closed=True),
        Contour("hole", [(5.004, 2.5)], [0], kind="point"),
        Contour("hole2", [(5, 20)], [0], kind="point"),
    ]
    fixed, report = fix_contours(contours)
    assert fixed[0].points[0][1] == fixed[1].points[0][1]  # basi collineari
    assert fixed[2].points[0][0] == fixed[3].points[0][0]  # centri allineati in X


def test_large_error_is_not_touched():
    # un segmento inclinato di 3° è una scelta progettuale, non un errore
    pts = [(0, 0), (20, 0), (20 + 10 * math.cos(math.radians(93)), 10 * math.sin(math.radians(93))),
           (0, 10)]
    fixed, _ = fix_contours([Contour("p", pts, [0] * 4, closed=True)])
    (x1, y1), (x2, y2) = fixed[0].points[1:3]
    assert abs(math.degrees(math.atan2(y2 - y1, x2 - x1)) - 93) < 0.01


def test_dxf_roundtrip(tmp_path, dirty_contours):
    doc = ezdxf.new()
    msp = doc.modelspace()
    for c in dirty_contours:
        msp.add_lwpolyline([(x, y, b) for (x, y), b in zip(c.points, c.bulges)],
                           format="xyb", close=True)
    # cava specchiata (estrusione -Z), come capita con i comandi SPECCHIO
    mirrored = msp.add_lwpolyline([(-x, y, -b) for (x, y), b in zip(dirty_contours[1].points,
                                                                   dirty_contours[1].bulges)],
                                  format="xyb", close=True, dxfattribs={"extrusion": (0, 0, -1)})
    mirrored.translate(100, 0, 0)
    arc = msp.add_arc((0, 60), 5, 0, 90.01)
    line = msp.add_line((5, 60), (5.002, 50))
    src = tmp_path / "matrice.dxf"
    doc.saveas(src)

    assert main([str(src), "-o", str(tmp_path / "out.dxf")]) == 0
    out = ezdxf.readfile(tmp_path / "out.dxf")
    contours, _ = dxf_io.extract(out.modelspace())
    for c in contours:
        if c.kind == "polyline":
            for dx, dy in segment_angles(c.points, c.bulges, True):
                angle = math.degrees(math.atan2(dy, dx)) % 180
                assert min(abs(angle - a) for a in (0, 30, 45, 90, 120, 135, 180)) < 1e-9
            for err, is_arc in tangent_errors(c.points, c.bulges):
                if is_arc:
                    assert err < 1e-9
    # linea resa verticale e arco reso tangente ad essa
    line_c = next(c for c in contours if c.kind == "line")
    assert line_c.points[0][0] == line_c.points[1][0]
    arc_e = next(e for e in out.modelspace() if e.dxftype() == "ARC")
    assert abs((arc_e.dxf.start_angle + 180) % 360 - 180) < 1e-9


def test_exact_design_angle_wins_over_mean():
    # lato lungo a 13° esatti (valore di progetto) e lato interno a 13.068°:
    # il gruppo deve prendere 13°, non la media pesata
    def seg(x, y, deg, length):
        r = math.radians(deg)
        return Contour(f"s{x}", [(x, y), (x + length * math.cos(r), y + length * math.sin(r))],
                       [0, 0], kind="line")
    fixed, report = fix_contours([seg(0, 0, 13.0, 36.6), seg(0, 5, 13.068474, 17.7)],
                                 Params(wall_step=0))
    for fc in fixed:
        (x0, y0), (x1, y1) = fc.points
        assert abs(math.degrees(math.atan2(y1 - y0, x1 - x0)) - 13.0) < 1e-9
    assert report.max_displacement < 0.02


def test_coincidental_arc_endpoint_alignment_is_ignored():
    # due raccordi diversi con estremi a X quasi uguale (0.6 micron) per caso:
    # non devono essere forzati in colonna né bloccare la correzione
    outer_pts, outer_b = fillet_polygon([(0, 0), (40, 0), (40, 20), (0, 20)], 3.0)
    inner_pts, inner_b = fillet_polygon([(1.0006, 3), (37, 3), (37, 17), (1.0006, 17)], 2.0)
    contours = [Contour("o", outer_pts, outer_b, closed=True),
                Contour("i", inner_pts, inner_b, closed=True)]
    fixed, report = fix_contours(contours, Params(wall_step=0))
    assert not report.warnings
    assert report.max_displacement < 1e-6


def test_exact_geometry_is_not_moved():
    # contorno esatto (spessore 1.2) + un lato interno storto di 0.05°:
    # si deve muovere solo il lato storto, lo spessore 1.2 resta identico
    r = math.radians(0.05)
    outer = Contour("o", [(0, 0), (40, 0), (40, 20), (0, 20)], [0] * 4, closed=True)
    inner = Contour("i", [(1.2, 1.2), (38.8, 1.2 + 37.6 * math.tan(r)), (38.8, 18.8), (1.2, 18.8)],
                    [0] * 4, closed=True)  # fondo interno inclinato
    fixed, report = fix_contours([outer, inner])
    pts = fixed[1].points
    assert fixed[0].points == outer.points  # contorno esatto intatto
    assert [p[0] for p in pts] == [1.2, 38.8, 38.8, 1.2]  # spessori laterali 1.2 intatti
    assert pts[2:] == [(38.8, 18.8), (1.2, 18.8)]  # lato superiore intatto
    assert pts[0][1] == pts[1][1]  # fondo reso orizzontale


def test_wall_thickness_snapped_to_one_decimal_and_radii_kept():
    # fondo inclinato a 13° spesso 1.6, con la linea interna disegnata storta
    # (spessore da 1.6 a 1.62) e raccordi R0.5 esterni / R0.3 interni:
    # lo spessore torna 1.6 costante e i raggi restano esattamente R0.5 / R0.3
    a = math.radians(13)
    t = 1.6
    outer = [(0, 40 * math.tan(a)), (40, 0), (40, 20), (0, 20)]
    inner_y = lambda x: (40 - x) * math.tan(a) + t / math.cos(a)
    inner = [(1.2, inner_y(1.2) + 0.02), (38.8, inner_y(38.8)), (38.8, 18.8), (1.2, 18.8)]
    o_pts, o_b = fillet_polygon(outer, 0.5)
    i_pts, i_b = fillet_polygon(inner, 0.3)
    fixed, report = fix_contours([Contour("o", o_pts, o_b, closed=True),
                                  Contour("i", i_pts, i_b, closed=True)])
    assert not report.warnings and not report.notes
    assert report.walls_changed >= 1
    o, i = fixed
    n = (math.sin(a), math.cos(a))  # normale al fondo (direzione -13°)
    off = lambda p: n[0] * p[0] + n[1] * p[1]
    for p in i.points[1:3]:  # il segmento 1 è il fondo in entrambi i contorni
        assert abs(off(p) - off(o.points[1]) - t) < 1e-9
    # il contorno esterno era esatto: non si muove
    assert all(math.dist(p, q) < 1e-12 for p, q in zip(o.points, o_pts))
    for fc, r_expected in ((o, 0.5), (i, 0.3)):
        pts, bs = fc.points, fc.bulges
        for k, b in enumerate(bs):
            if b:
                (x0, y0), (x1, y1) = pts[k], pts[(k + 1) % len(pts)]
                r = math.hypot(x1 - x0, y1 - y0) / (2 * abs(math.sin(2 * math.atan(b))))
                assert abs(r - r_expected) < 1e-9
