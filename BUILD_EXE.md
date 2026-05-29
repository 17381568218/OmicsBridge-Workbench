# OmicsBridge Metabolomics Desktop Build

Desktop entry:

```powershell
python -B metabolomics_gui.py
```

First-stage analysis input:

- CSV matrix with one feature ID column, for example `Compound_ID`
- Sample columns such as `Raw_20250528-HXY-NEG-C1-1`

Outputs:

- `differential_<GroupA>_vs_<GroupB>.csv`
- `sample_design.csv`
- `qc_summary.csv`
- `analysis_report.txt`
- `pca.png`
- `volcano.png`
- `heatmap.png`

PyInstaller build command:

```powershell
pyinstaller --noconfirm --windowed --name OmicsBridgeMetabolomics metabolomics_gui.py
```

If pandas, numpy, or matplotlib hidden imports are missed by the local Python
environment, use:

```powershell
pyinstaller --noconfirm --windowed --name OmicsBridgeMetabolomics `
  --hidden-import pandas --hidden-import numpy --hidden-import matplotlib `
  metabolomics_gui.py
```

Recommended packaging direction:

- Keep `metabolomics.py` as the reusable analysis engine.
- Keep `metabolomics_gui.py` as the desktop GUI entry.
- Add proteomics and transcriptomics as separate engine modules later.
- Only keep the Flask web app as an optional developer/debug interface.
