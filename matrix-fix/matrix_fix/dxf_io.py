"""Lettura/scrittura delle entità DXF/DWG verso il modello di ``core``.

Supporta LWPOLYLINE, POLYLINE 2D, LINE, ARC e CIRCLE. Le coordinate vengono
portate in WCS (gestendo anche le entità specchiate con estrusione 0,0,-1) e
riscritte nel loro OCS originale. I file DWG richiedono l'ODA File Converter
(gratuito) installato: ezdxf lo usa per convertire DWG <-> DXF.
"""
from __future__ import annotations

import math
from pathlib import Path

import ezdxf
from ezdxf.math import Vec3

from .core import Contour, FixedContour

SUPPORTED = {"LWPOLYLINE", "POLYLINE", "LINE", "ARC", "CIRCLE"}


def load(path: str | Path):
    path = Path(path)
    if path.suffix.lower() == ".dwg":
        from ezdxf.addons import odafc
        return odafc.readfile(str(path))
    return ezdxf.readfile(str(path))


def save(doc, path: str | Path) -> None:
    path = Path(path)
    if path.suffix.lower() == ".dwg":
        from ezdxf.addons import odafc
        odafc.export_dwg(doc, str(path), replace=True)
    else:
        doc.saveas(str(path))


def _flip(entity) -> bool | None:
    """True se l'OCS è specchiato (estrusione -Z), None se non è nel piano XY."""
    ext = Vec3(entity.dxf.get("extrusion", (0, 0, 1))).normalize()
    if ext.isclose(Vec3(0, 0, 1)):
        return False
    if ext.isclose(Vec3(0, 0, -1)):
        return True
    return None


def _to_wcs(entity, x: float, y: float) -> tuple[float, float]:
    p = entity.ocs().to_wcs(Vec3(x, y, 0))
    return (p.x, p.y)


def _to_ocs(entity, x: float, y: float) -> tuple[float, float]:
    p = entity.ocs().from_wcs(Vec3(x, y, 0))
    return (p.x, p.y)


def extract(layout, layers: set[str] | None = None) -> tuple[list[Contour], dict[str, int]]:
    """Estrae i contorni dal layout; restituisce anche il conteggio delle entità saltate."""
    contours: list[Contour] = []
    skipped: dict[str, int] = {}

    def skip(reason: str) -> None:
        skipped[reason] = skipped.get(reason, 0) + 1

    for e in layout:
        kind = e.dxftype()
        if layers is not None and e.dxf.get("layer", "0") not in layers:
            continue
        if kind not in SUPPORTED:
            skip(kind)
            continue
        if kind == "LINE":
            s, t = e.dxf.start, e.dxf.end
            contours.append(Contour(e, [(s.x, s.y), (t.x, t.y)], [0.0, 0.0], kind="line"))
            continue
        flip = _flip(e)
        if flip is None:
            skip(f"{kind} fuori dal piano XY")
            continue
        sign = -1.0 if flip else 1.0
        if kind == "LWPOLYLINE":
            pts = e.get_points("xyb")
            contours.append(Contour(
                e, [_to_wcs(e, x, y) for x, y, _ in pts], [sign * b for _, _, b in pts],
                closed=e.closed))
        elif kind == "POLYLINE":
            if not e.is_2d_polyline or e.dxf.flags & (2 | 4):
                skip("POLYLINE 3D/mesh/smussata")
                continue
            verts = list(e.vertices)
            contours.append(Contour(
                e, [_to_wcs(e, v.dxf.location.x, v.dxf.location.y) for v in verts],
                [sign * v.dxf.get("bulge", 0.0) for v in verts], closed=e.is_closed))
        elif kind == "ARC":
            c, r = e.dxf.center, e.dxf.radius
            sa, ea = math.radians(e.dxf.start_angle), math.radians(e.dxf.end_angle)
            sweep = (ea - sa) % (2 * math.pi)
            if sweep < 1e-12:
                skip("ARC di ampiezza nulla")
                continue
            p = _to_wcs(e, c.x + r * math.cos(sa), c.y + r * math.sin(sa))
            q = _to_wcs(e, c.x + r * math.cos(ea), c.y + r * math.sin(ea))
            contours.append(Contour(e, [p, q], [sign * math.tan(sweep / 4), 0.0], kind="arc"))
        elif kind == "CIRCLE":
            c = e.dxf.center
            contours.append(Contour(e, [_to_wcs(e, c.x, c.y)], [0.0], kind="point"))
    return contours, skipped


def apply(fixed: list[FixedContour], contours: list[Contour]) -> int:
    """Scrive le correzioni nelle entità; restituisce il numero di entità eliminate."""
    deleted = 0
    for fc, c in zip(fixed, contours):
        e = c.key
        kind = e.dxftype()
        if fc.deleted:  # entità degenere (lunghezza nulla entro la tolleranza)
            e.destroy()
            deleted += 1
            continue
        if kind == "LINE":
            (x0, y0), (x1, y1) = fc.points
            e.dxf.start = Vec3(x0, y0, e.dxf.start.z)
            e.dxf.end = Vec3(x1, y1, e.dxf.end.z)
            continue
        sign = -1.0 if _flip(e) else 1.0
        if kind == "LWPOLYLINE":
            old = e.get_points("xyseb")
            pts = []
            for (x, y), b, src in zip(fc.points, fc.bulges, fc.src_index):
                ox, oy = _to_ocs(e, x, y)
                pts.append((ox, oy, old[src][2], old[src][3], sign * b))
            e.set_points(pts, format="xyseb")
            e.closed = fc.closed
        elif kind == "POLYLINE":
            verts = list(e.vertices)
            keep = set(fc.src_index)
            for i in range(len(verts) - 1, -1, -1):
                if i not in keep:
                    e.delete_vertices(i, 1)
            for v, (x, y), b in zip(e.vertices, fc.points, fc.bulges):
                ox, oy = _to_ocs(e, x, y)
                v.dxf.location = Vec3(ox, oy, v.dxf.location.z)
                v.dxf.bulge = sign * b
            e.close(fc.closed)
        elif kind == "ARC":
            (px, py), (qx, qy) = (_to_ocs(e, *pt) for pt in fc.points)
            b = sign * fc.bulges[0]
            if b < 0:  # ARC è sempre antiorario nel suo OCS
                (px, py), (qx, qy), b = (qx, qy), (px, py), -b
            chord = math.hypot(qx - px, qy - py)
            ux, uy = (qx - px) / chord, (qy - py) / chord
            h = chord * (1 - b * b) / (4 * b)  # distanza (con segno) centro-corda
            cx, cy = (px + qx) / 2 - uy * h, (py + qy) / 2 + ux * h
            e.dxf.center = Vec3(cx, cy, e.dxf.center.z)
            e.dxf.radius = math.hypot(px - cx, py - cy)
            e.dxf.start_angle = math.degrees(math.atan2(py - cy, px - cx))
            e.dxf.end_angle = math.degrees(math.atan2(qy - cy, qx - cx))
        elif kind == "CIRCLE":
            cx, cy = _to_ocs(e, *fc.points[0])
            e.dxf.center = Vec3(cx, cy, e.dxf.center.z)
    return deleted
