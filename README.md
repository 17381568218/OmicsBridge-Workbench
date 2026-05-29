# OmicsBridge Workbench

A comprehensive web-based multi-omics annotation and analysis platform integrating CCS (Collision Cross Section) annotation, HMDB source classification, metabolomics statistical analysis, and compound identification post-processing.

## Features

### 1. CCS Annotation Pipeline (Tab 1)
- Dual ionization mode support (POS + NEG)
- Automatic CCS database matching (650K+ entries from AllCCS2, CCSbase, PNNL, NORMAN SLE)
- Adduct determination via m/z and neutral mass
- HMDB biological source annotation (Origin, Biospecimen, Tissue)
- Chemical class annotation from HMDB taxonomy (145K+ compounds)
- Compound-level deduplication and HMDB-level deduplication
- Score / Fragmentation Score / dCCS multi-criteria filtering
- Origin and Tissue priority selection
- COM raw abundance matrix extraction
- POS+NEG merged output with Compound ID deduplication

### 2. Metabolomics Analysis (Tab 2)
- Feature-by-sample matrix upload (CSV)
- Automatic sample group detection
- QC-RSD quality assessment
- PCA (Principal Component Analysis) with interactive SVG plots
- Differential analysis (Welch's t-test + permutation p-values)
- Volcano plot visualization
- Top feature heatmap
- Exportable result tables

### 3. HMDB Source Query (Tab 3)
- Query Origin, Biospecimen, Tissue, Cellular location by HMDB ID
- Paste text or import CSV/XLSX files

### 4. CCS Value Query (Tab 4)
- Look up experimental and predicted CCS values by HMDB ID + Adduct
- Paste text or import CSV/XLSX files

### 5. Origin Query (Tab 5)
- Batch query biological origin (Endogenous, Food, Drug, Environmental, etc.) by HMDB ID
- Paste text or import CSV/XLSX files

### 6. Chemical Class Query (Tab 6)
- Batch query chemical taxonomy class by HMDB ID
- Paste text or import CSV/XLSX files

## System Requirements

- **Python** 3.10+ 
- **Operating System**: Windows / macOS / Linux
- **RAM**: 4GB minimum (8GB recommended for large datasets)
- **Disk**: ~300MB (including 38MB CCS database + 5MB chemical class data)

## Installation

```bash
# Clone repository
git clone https://github.com/YOUR_USERNAME/omicsbridge-workbench.git
cd omicsbridge-workbench

# Install dependencies
pip install -r requirements.txt
```

## Quick Start

```bash
# Start the web server
python run.py
```

Open your browser to `http://localhost:5000`

### Public Access (Optional)

To share with collaborators over the internet:

```bash
# Using Cloudflare Tunnel (free, temporary URL)
cloudflared tunnel --url http://localhost:5000 --protocol http2
```

See `start_tunnel.bat` for a packaged startup script.

## Usage Guide

### CCS Annotation Pipeline Workflow

1. **Upload Files**: Navigate to Tab 1, select your POS COM and IDEN CSV files. Optionally add NEG mode files.
2. **Configure Filters**: Set Score/Fragmentation Score/dCCS cutoffs. Select Origins to keep and Tissues to prioritize.
3. **Run Pipeline**: Click "Run". The pipeline will:
   - Load the CCS database (650K+ compounds)
   - Annotate each compound with CCS values (Experimental > Predicted priority)
   - Add HMDB biological source and chemical class
   - Filter by Score, Fragmentation Score, dCCS
   - Deduplicate by Compound ID with loser retry
   - Merge COM raw abundance data
4. **Download Results**: 
   - `output_POS.csv` / `output_NEG.csv` — Single mode results
   - `output_merged.csv` — Combined POS+NEG with Compound ID dedup
   - `output_merged_raw_matrix.csv` — Compound ID × Raw abundance matrix

### Metabolomics Analysis Workflow

1. **Upload Matrix**: Navigate to Tab 2, select a CSV where rows are features and columns are samples.
2. **Set Groups** (optional): Define Group A and Group B labels for differential analysis. Auto-detected if left blank.
3. **Run Analysis**: View PCA scores plot, volcano plot, differential features table, heatmap, and QC metrics.
4. **Download**: Export result tables and figures.

## Input File Formats

### COM File (Compound Measurement)
Multi-row header CSV exported from Progenesis QI:
- Row 1: Section labels ("Normalised abundance", "Raw abundance")
- Row 2: Sample group labels (C, H, QC)
- Row 3: Sample names
- Column 0: Compound name
- Column 5: CCS value
- Raw abundance: columns after "Raw abundance" label

### IDEN File (Compound Identification)
Single or two-row header CSV:
- Column 0: Compound name
- Column 1: Compound ID (HMDB or Metlin)
- Column 5: Score
- Column 6: Fragmentation Score

### Metabolomics Matrix
Standard CSV format:
- First column: Feature/Compound ID
- Remaining columns: Sample abundance values
- Sample names should contain group identifiers (e.g., "C1-1", "H2-1", "QC-1")

## Architecture

```
omicsbridge_workbench/
├── app.py                 # Flask web server + REST API
├── processor.py           # CCS annotation core logic
├── metabolomics.py        # Metabolomics analysis engine
├── run.py                 # Production launcher (waitress WSGI)
├── requirements.txt       # Python dependencies
├── data/
│   ├── ccs_data.db.gz     # CCS database (650K+ entries)
│   └── hmdb_chemclass.json # Chemical classification data
├── static/
│   ├── app.js             # Frontend logic
│   └── style.css          # UI styles
└── templates/
    └── index.html         # Single-page web application
```

### API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/run` | POST | Launch CCS annotation pipeline |
| `/api/metabolomics/run` | POST | Launch metabolomics analysis |
| `/api/status/<task_id>` | GET | Poll task status and logs |
| `/api/download/<task_id>/<filename>` | GET | Download output files |
| `/api/hmdb_query` | POST | HMDB source annotation query |
| `/api/ccs_query` | POST | CCS value lookup |
| `/api/origin_query` | POST | Biological origin query |
| `/api/chemclass_query` | POST | Chemical class query |
| `/api/origins_tissues` | GET | Get available filter options |

## Databases

The CCS database integrates data from:
- **AllCCS2** (126K HMDB-annotated compounds)
- **CCSbase** (experimental CCS values)
- **PNNL** (high-quality experimental CCS)
- **NORMAN SLE** (environmental compounds)

The HMDB biological source index covers 298K+ metabolites with:
- Biological origin classification (Endogenous, Food, Drug, Environmental, Microbial, Synthetic)
- Biospecimen and tissue locations
- Chemical taxonomy class (from HMDB 5.0 XML, 145K+ compounds)

## Desktop Version

A standalone desktop GUI (`metabolomics_gui.py`) is also available for offline use. Build with:

```bash
pyinstaller --onedir --noconsole metabolomics_gui.py
```

See `BUILD_EXE.md` for detailed build instructions.

## License

This project is for academic and research use.

## Citation

If you use this software in your research, please cite the underlying databases:
- AllCCS2: [allccs.zhulab.cn](http://allccs.zhulab.cn)
- HMDB: [hmdb.ca](https://hmdb.ca)
- CCSbase: [ccsbase.net](http://ccsbase.net)
- PNNL CCS Database
