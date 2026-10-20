@echo off
rem Trascina un file .dwg o .dxf su questo file per correggerlo.
rem Il risultato viene scritto accanto all originale come <nome>_fix.<ext>
cd /d "%~dp0"
if "%~1"=="" (
  echo Trascina un file .dwg o .dxf su matrix-fix.bat
) else (
  python -m matrix_fix "%~1" -v %2 %3 %4 %5 %6 %7 %8 %9
)
pause
