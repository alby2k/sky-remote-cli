"""Riga di comando: python -m matrix_fix ingresso.dxf [-o uscita.dxf] [opzioni]."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import dxf_io
from .core import Params, fix_contours


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="matrix_fix",
        description="Corregge micro-errori geometrici (ortogonalità, parallelismo, "
                    "allineamenti, tangenze, aperture) nei profili DXF/DWG di matrici.")
    p.add_argument("input", type=Path, help="file .dxf o .dwg")
    p.add_argument("-o", "--output", type=Path,
                   help="file di uscita (default: <nome>_fix.<ext>)")
    p.add_argument("--lin-tol", type=float, default=0.02,
                   help="tolleranza lineare in mm (default 0.02)")
    p.add_argument("--angle-tol", type=float, default=0.1,
                   help="errore angolare sempre corretto, in gradi (default 0.1)")
    p.add_argument("--max-angle", type=float, default=1.0,
                   help="errore angolare massimo correggibile sui segmenti corti (default 1.0)")
    p.add_argument("--angle-step", type=float, default=15.0,
                   help="griglia angolare in gradi; 0 = solo 0/90 (default 15)")
    p.add_argument("--tangent-tol", type=float, default=1.0,
                   help="scarto oltre il quale una giunzione arco/segmento è uno spigolo (default 1.0)")
    p.add_argument("--merge-tol", type=float,
                   help="distanza sotto la quale due vertici vengono fusi (default = lin-tol)")
    p.add_argument("--max-move", type=float,
                   help="spostamento massimo ammesso di un vertice in mm (default 5 x lin-tol)")
    p.add_argument("--no-align-vertices", action="store_true",
                   help="non allineare in X/Y i vertici quasi allineati")
    p.add_argument("--layers", help="elabora solo questi layer (separati da virgola)")
    p.add_argument("--include-blocks", action="store_true",
                   help="corregge anche le definizioni di blocco (ognuna indipendente)")
    p.add_argument("--check", action="store_true",
                   help="analizza e riporta gli errori senza scrivere il file")
    p.add_argument("--force", action="store_true",
                   help="scrive anche se qualche vertice supera --max-move")
    p.add_argument("-v", "--verbose", action="store_true", help="elenca ogni correzione")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    params = Params(lin_tol=args.lin_tol, angle_tol=args.angle_tol, max_angle=args.max_angle,
                    angle_step=args.angle_step, tangent_tol=args.tangent_tol,
                    merge_tol=args.merge_tol, max_move=args.max_move,
                    align_vertices=not args.no_align_vertices)
    layers = set(args.layers.split(",")) if args.layers else None

    try:
        doc = dxf_io.load(args.input)
    except Exception as exc:  # file illeggibile o ODA File Converter mancante
        print(f"Impossibile leggere {args.input}: {exc}", file=sys.stderr)
        return 2

    layouts = [("Model", doc.modelspace())]
    if args.include_blocks:
        layouts += [(f"blocco {b.name}", b) for b in doc.blocks
                    if not b.is_any_layout and not b.name.startswith("*")]

    blocked = False
    for name, layout in layouts:
        contours, skipped = dxf_io.extract(layout, layers)
        if not contours:
            continue
        fixed, report = fix_contours(contours, params)
        print(f"== {name}: {len(contours)} entità")
        print(report.summary())
        for reason, count in skipped.items():
            print(f"  ignorate: {count} x {reason}")
        if args.verbose:
            for line in report.details:
                print("  -", line)
        for note in report.notes:
            print("NOTA:", note)
        for warning in report.warnings:
            print("ATTENZIONE:", warning)
        if report.warnings and not args.force:
            blocked = True
            continue
        if not args.check:
            deleted = dxf_io.apply(fixed, contours)
            if deleted:
                print(f"  entità degeneri eliminate: {deleted}")

    if args.check:
        return 0
    if blocked:
        print("File NON scritto: rivedere le tolleranze o usare --force.", file=sys.stderr)
        return 1
    out = args.output or args.input.with_name(f"{args.input.stem}_fix{args.input.suffix}")
    dxf_io.save(doc, out)
    print(f"Scritto {out}")
    return 0
