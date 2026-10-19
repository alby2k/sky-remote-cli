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
    assert abs(arc_e.dxf.start_angle % 360) < 1e-9
