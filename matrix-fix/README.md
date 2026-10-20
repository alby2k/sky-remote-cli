# matrix-fix

Correzione automatica dei micro-errori geometrici nei disegni DXF/DWG delle
matrici di estrusione (profili di facciata), prima di importarli come schizzo
in Autodesk Inventor.

Errori corretti (tutti nell'ordine dei centesimi/millesimi di mm o di grado):

| Errore | Cosa fa l'algoritmo |
|---|---|
| Lati quasi a 90° | gruppi di direzioni quasi ortogonali resi esattamente ortogonali |
| Linee di spessore quasi parallele | direzioni quasi uguali (anche su polilinee diverse) unificate |
| Linee quasi orizzontali/verticali | agganciate a 0°/90° in modo *esatto* (Y o X identiche bit per bit) |
| Angoli quasi "tondi" (45°, 30°, 60°…) | agganciati alla griglia angolare (`--angle-step`, default 15°) |
| Lati quasi allineati | segmenti paralleli con offset quasi uguale resi collineari |
| Vertici quasi allineati | vertici/centri di fori con X (o Y) quasi uguale allineati |
| Raccordi non perfettamente tangenti | angolo dell'arco e corda ricalcolati per tangenza esatta, raggio di progetto conservato |
| Spessori cartella "storti" (es. 1,567 mm) | linee parallele a ≤ 4 mm portate a distanza multipla di 0,1 mm, costante |
| Quote orizzontali/verticali non al decimo (es. 7,954) | tutte le distanze fra linee orizzontali (e fra verticali) portate al decimo |
| Contorni non chiusi, micro-aperture | vertici più vicini di `--merge-tol` fusi, polilinee chiuse |
| Micro-segmenti, micro-archi | eliminati / raddrizzati |

## Come funziona

1. **Topologia** – tutti i vertici (di polilinee, linee, archi) a distanza minore
   di `merge-tol` diventano un unico nodo condiviso: si chiudono le aperture fra
   entità e spariscono i segmenti di lunghezza quasi nulla.
2. **Direzioni** – i segmenti, dal più lungo al più corto, vengono raggruppati per
   direzione. Un segmento entra in un gruppo se lo scarto angolare è ≤ `angle-tol`
   **oppure** se lo spostamento laterale che produce (`lunghezza × sin(scarto)`) è
   ≤ `lin-tol` (con un tetto di `max-angle`). Questo secondo criterio è
   importante: un errore di 0,01 mm su un segmento corto di 0,5 mm vale più di 1°,
   mentre su un lato di 50 mm vale 0,01°. Ogni gruppo prende l'angolo "di
   progetto" se una sua linea ce l'ha già esatto (es. 13,000° o 29,000°),
   altrimenti la direzione media pesata sulle lunghezze; poi viene agganciato a
   0°/90°/griglia se vicino (un angolo tondo esatto come 29° non viene mai
   portato a 30°), e i gruppi restanti vengono resi ortogonali a quelli già
   fissati se ci sono vicini.
3. **Tangenze** – in ogni giunzione arco/segmento (o arco/arco) quasi tangente
   (≤ `tangent-tol`), si percorre la catena di entità tangenti: fra due segmenti
   di direzione ormai esatta la somma degli angoli al centro degli archi deve
   essere uguale alla rotazione fra le due direzioni; la piccola differenza viene
   ripartita sugli archi e da qui si ricava la direzione esatta di ogni corda.
   L'allineamento dei vertici riguarda solo spigoli fra segmenti e centri dei
   fori: gli estremi dei raccordi sono esclusi (li posiziona la tangenza).
   I raggi dei raccordi restano quelli di progetto (arrotondati al centesimo se
   ci sono vicini): con direzione della corda nota, anche il raggio è un
   vincolo lineare.
4. **Distanze al decimo** (`--wall-step`, default 0,1 mm):
   - **linee orizzontali e verticali**: tutte le distanze fra loro sono multipli di
     0,1 mm. Le linee di ciascuna direzione vengono messe su una griglia di passo
     0,1 la cui origine è quella su cui cade già la maggior parte delle linee
     (pesate per lunghezza): le linee già in passo restano *ferme*, si spostano solo
     quelle fuori passo. Una linea esattamente a metà fra due decimi (scarto
     ≈ 0,05) non si può decidere: resta com'è ed è segnalata come NOTA;
   - **spessori cartella**: due segmenti paralleli di qualunque direzione,
     affacciati, a non più di `--wall-max` (default 4 mm) e senza altri segmenti
     fra loro, vengono portati a distanza multipla di 0,1 mm.
5. **Vincoli lineari** – fissate le direzioni, ogni condizione diventa un'equazione
   *lineare* nelle coordinate dei vertici, ad esempio per un segmento A→B di
   direzione con normale `n`: `n · (B − A) = 0`; per due segmenti collineari:
   `n · (A₂ − A₁) = 0`; per due vertici allineati in verticale: `x₁ − x₂ = 0`.
6. **Minimo spostamento** – si risolve il sistema cercando lo spostamento di
   norma minima dei vertici (minimi quadrati vincolati, `scipy lsqr`): il disegno
   si muove solo di quanto serve, tipicamente pochi micron. I vertici fra linee
   già esatte (e fra raccordi rimasti invariati) pesano molto di più, così la
   correzione ricade sulle linee davvero sbagliate e le quote già corrette
   (es. spessori di 1,200 mm) non cambiano; ciò che si muoverebbe meno di
   1e-9 mm resta identico all'originale. Le uguaglianze su X/Y
   vengono infine imposte bit per bit, così Inventor vede linee *esattamente*
   orizzontali/verticali.
7. **Sicurezza** – se un vertice dovesse spostarsi più di `max-move`
   (default 5 × `lin-tol`) il file non viene scritto e viene segnalato il punto:
   vuol dire che un vincolo non è un micro-errore ma una scelta progettuale.
   Gli scarti oltre `max-angle` (es. un lato a 93°) non vengono toccati.
   Se un vincolo facoltativo (spessore, raggio, collinearità, allineamento)
   richiederebbe spostamenti eccessivi o contraddice gli altri, viene scartato e
   segnalato come NOTA (priorità: distanze al decimo, collinearità, raggi,
   allineamenti). Se un raggio viene lasciato libero solo per far chiudere una
   catena di raccordi (scarti di milionesimi di mm) è riportato nel dettaglio.

## Uso su Windows (passo passo)

1. **Python**: installa Python 3.11 o successivo da [python.org](https://www.python.org/downloads/)
   spuntando *"Add python.exe to PATH"* nella prima schermata dell'installer.
2. **ODA File Converter** (per leggere/scrivere i DWG): scaricalo gratis da
   [opendesign.com](https://www.opendesign.com/guestfiles/oda_file_converter) e installalo
   normalmente. Il programma lo trova da solo in `C:\Program Files\ODA\...`;
   altrimenti si indica con `--oda "C:\percorso\ODAFileConverter.exe"`.
3. **Codice**: scarica la cartella `matrix-fix` (da GitHub: branch
   `claude/dwg-geometry-error-correction-7y1w4m` → *Code* → *Download ZIP*, oppure
   `git clone -b claude/dwg-geometry-error-correction-7y1w4m https://github.com/alby2k/sky-remote-cli.git`).
4. **Librerie**: apri il *Prompt dei comandi* nella cartella `matrix-fix` e lancia
   una volta sola:
   ```
   python -m pip install -r requirements.txt
   ```
5. **Uso**: trascina il file `.dwg` sopra `matrix-fix.bat`, oppure dal prompt:
   ```
   python -m matrix_fix "C:\cartella\Matrice.dwg" --check -v   (solo analisi)
   python -m matrix_fix "C:\cartella\Matrice.dwg" -v           (scrive Matrice_fix.dwg)
   ```
   Il file originale non viene mai modificato.

## Installazione

```
pip install -r requirements.txt
```

In alternativa all'ODA File Converter, i DWG si possono convertire in DXF anche con
`dwg2dxf` di [LibreDWG](https://github.com/LibreDWG/libredwg) (open source).

Per leggere/scrivere direttamente i **DWG** serve l'[ODA File Converter](https://www.opendesign.com/guestfiles/oda_file_converter)
(gratuito); in alternativa salvare il DWG come DXF da AutoCAD/DWG TrueView.

## Uso

```
python -m matrix_fix matrice.dwg                 # scrive matrice_fix.dwg
python -m matrix_fix matrice.dxf -o pulita.dxf -v # elenca ogni correzione
python -m matrix_fix matrice.dxf --check          # solo analisi, nessun file scritto
```

Opzioni principali:

| Opzione | Default | Significato |
|---|---|---|
| `--lin-tol` | 0.02 mm | errore lineare massimo considerato micro-errore |
| `--angle-tol` | 0.1° | errore angolare sempre corretto |
| `--max-angle` | 0.25° | errore angolare massimo correggibile (segmenti corti) |
| `--angle-step` | 15° | griglia angolare; `0` = solo 0°/90° |
| `--tangent-tol` | 1.0° | oltre questo scarto una giunzione è uno spigolo vivo |
| `--merge-tol` | = lin-tol | distanza di fusione dei vertici |
| `--max-move` | 5 × lin-tol | spostamento massimo ammesso di un vertice |
| `--wall-step` | 0.1 mm | passo di distanze e spessori; `0` = non toccarli |
| `--wall-tol` | 0.049 mm | scarto massimo di una distanza dal valore arrotondato |
| `--wall-max` | 4 mm | linee parallele fino a questa distanza sono spessori cartella |
| `--no-align-vertices` | | non allineare i vertici in X/Y |
| `--layers A,B` | tutti | elabora solo alcuni layer |
| `--include-blocks` | | corregge anche le definizioni di blocco |

Entità gestite: LWPOLYLINE, POLYLINE 2D, LINE, ARC, CIRCLE (anche specchiate,
estrusione 0,0,−1). SPLINE ed ELLISSI vengono ignorate e segnalate.

## Test

```
pip install pytest
python -m pytest -q
```

I test costruiscono un profilo con smusso a 45°, cava interna con raccordi,
parete a 30° e contorni fatti di LINE separate, lo "sporcano" con errori casuali
di ±5 µm e verificano che dopo la correzione direzioni, parallelismi, ortogonalità
e tangenze siano esatti (scarto < 1e-9°) con spostamenti sotto 0,02 mm.
