"""
Core CCS annotation processing logic (extracted from ccs_annotate_gui.py).
No tkinter dependencies — usable from Flask backend.
"""
import csv, re, json, sqlite3, gzip, shutil, tempfile, os, sys, time
from collections import defaultdict

# ============================================================
# Path helpers
# ============================================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def get_db_path():
    for p in [os.path.join(BASE_DIR, 'data', 'ccs_data.db.gz'),
              os.path.join(BASE_DIR, 'data', 'ccs_data.db')]:
        if os.path.exists(p): return p
    return ''

def get_data_path(filename):
    for d in [os.path.join(BASE_DIR, 'data')]:
        p = os.path.join(d, filename)
        if os.path.exists(p): return p
    return ''

def open_db():
    db_path = get_db_path()
    if not db_path: raise FileNotFoundError("CCS database not found")
    db_file = db_path
    if db_path.endswith('.gz'):
        td = os.path.join(tempfile.gettempdir(), 'ccs_web.db')
        if not os.path.exists(td) or os.path.getmtime(db_path) > os.path.getmtime(td):
            with gzip.open(db_path, 'rb') as gz:
                with open(td, 'wb') as out:
                    shutil.copyfileobj(gz, out)
        db_file = td
    conn = sqlite3.connect(db_file)
    conn.execute('PRAGMA query_only=ON')
    return conn

def load_chemclass():
    chemclass = {}
    json_path = get_data_path('hmdb_chemclass.json')
    if json_path and os.path.exists(json_path):
        with open(json_path, encoding='utf-8') as f:
            data = json.load(f)
            for hid, cc in data.items():
                chemclass[hid.upper()] = cc
    return chemclass

def load_kegg():
    kegg = {}
    json_path = get_data_path('hmdb_kegg.json')
    if json_path and os.path.exists(json_path):
        with open(json_path, encoding='utf-8') as f:
            data = json.load(f)
            for hid, kid in data.items():
                kegg[hid.upper()] = kid
    return kegg

# ============================================================
# Adduct / naming helpers
# ============================================================
def simplify_name(n):
    if not n: return ''
    return n.upper().replace(' ', '').replace('-', '').replace(',', '').replace('/', '').replace('(', '').replace(')', '').replace("'", '')

def norm_adduct(a):
    if not a: return ''
    return a.strip().replace('[','').replace(']','').replace('(','').replace(')','').upper()

ADDUCT_MASSES = {
    'M+H': 1.007276, 'M+NA': 22.989218, 'M+K': 38.963158,
    'M+NH4': 18.033823, 'M-H2O+H': -17.00274, 'M+NA-2H': 20.974666,
    'M+HCOO': 44.99765, 'M-H': -1.007276, 'M+CH3COO': 59.01385,
    'M+CL': 34.968403, '2M+H': 1.007276, '2M+NA': 22.989218, '2M-H': -1.007276,
}

def parse_adduct_mass_diff(a):
    a = a.upper().strip().replace('[', '').replace(']', '')
    if a in ADDUCT_MASSES: return ADDUCT_MASSES[a]
    for k, v in ADDUCT_MASSES.items():
        if a == k or a.replace('+', '+') == k: return v
    return None

def determine_best_adduct(neutral_mass, mz, adduct_list):
    if not neutral_mass or not mz: return None, None
    try: nm, mz_val = float(neutral_mass), float(mz)
    except: return None, None
    ba, be = None, float('inf')
    for ra in adduct_list:
        md = parse_adduct_mass_diff(ra.strip())
        if md is None: continue
        em = nm + md; ep = abs(mz_val - em) / em * 1e6
        if ep < be: be, ba = ep, ra.strip()
    return (ba, be) if ba and be < 50 else (None, None)

# ============================================================
# Ranking & Dedup functions
# ============================================================
def rank_compound_rows(rows, col_map, score_cutoff, frag_cutoff,
                       allowed_origins, priority_tissues, dccs_cutoff=None):
    score_col = col_map.get('score')
    frag_col = col_map.get('frag')
    dccs_col = col_map.get('dccs')
    origin_col = col_map.get('origin')
    bio_col = col_map.get('biospecimen')
    tissue_col = col_map.get('tissue')

    def _score(r):
        try: return float(r[score_col]) if score_col is not None else 0
        except: return 0
    def _frag(r):
        try: return float(r[frag_col]) if frag_col is not None else 0
        except: return 0
    def _dccs(r):
        try: return float(r[dccs_col]) if dccs_col is not None else float('inf')
        except: return float('inf')
    def _has_tissue(r):
        if not priority_tissues: return True
        bio_val = (r[bio_col] if bio_col is not None else '') or ''
        tis_val = (r[tissue_col] if tissue_col is not None else '') or ''
        combined = f"{bio_val}; {tis_val}".lower()
        return any(t.lower() in combined for t in priority_tissues)

    if score_cutoff is not None and score_col is not None:
        passed = [r for r in rows if _score(r) >= score_cutoff]
        if not passed: return None
    else:
        passed = list(rows)

    if frag_cutoff is not None and frag_col is not None:
        passed = [r for r in passed if _frag(r) >= frag_cutoff]
        if not passed: return None

    if allowed_origins and origin_col is not None:
        origin_passed = [r for r in passed
                         if any(a.strip().lower() == (r[origin_col] or '').strip().lower()
                                for a in allowed_origins)]
        if origin_passed:
            passed = origin_passed

    if dccs_cutoff is not None and dccs_col is not None:
        passed = [r for r in passed if _dccs(r) <= dccs_cutoff]
        if not passed: return None

    def sort_key(r):
        tissue_bonus = 0 if _has_tissue(r) or not priority_tissues else 1
        d = _dccs(r)
        if d < 5: return (tissue_bonus, 0, d, -_score(r))
        return (tissue_bonus, 1, -_score(r), d)

    passed.sort(key=sort_key)
    return passed

def compare_rows(row_a, row_b, col_map):
    dccs_col = col_map.get('dccs')
    score_col = col_map.get('score')
    def _dccs(r):
        try: return float(r[dccs_col]) if dccs_col is not None and dccs_col < len(r) else float('inf')
        except: return float('inf')
    def _score(r):
        try: return float(r[score_col]) if score_col is not None and score_col < len(r) else 0
        except: return 0
    da, db = _dccs(row_a), _dccs(row_b)
    a_good, b_good = da < 5, db < 5
    if a_good and not b_good: return row_a
    if b_good and not a_good: return row_b
    if a_good and b_good: return row_a if da <= db else row_b
    return row_a if _score(row_a) >= _score(row_b) else row_b

def hmdb_dedup(compound_candidates, col_map):
    id_col = 1
    assigned_ids = {}
    final_kept = {}
    retry_queue = []

    for compound, candidates in compound_candidates.items():
        best = candidates[0]
        cid = (best[id_col] or '').strip().upper() if id_col < len(best) else ''
        if not cid:
            final_kept[compound] = best
        elif cid not in assigned_ids:
            assigned_ids[cid] = (compound, best)
        else:
            existing_cpd, existing_row = assigned_ids[cid]
            winner = compare_rows(existing_row, best, col_map)
            if winner is existing_row:
                retry_queue.append((compound, 1))
            else:
                assigned_ids[cid] = (compound, best)
                retry_queue.append((existing_cpd, 1))

    retry_saved = 0
    retry_dropped = 0
    for compound, start_idx in retry_queue:
        candidates = compound_candidates[compound]
        found = False
        for idx in range(start_idx, len(candidates)):
            row = candidates[idx]
            cid = (row[id_col] or '').strip().upper() if id_col < len(row) else ''
            if not cid:
                final_kept[compound] = row
                retry_saved += 1
                found = True; break
            elif cid not in assigned_ids:
                assigned_ids[cid] = (compound, row)
                final_kept[compound] = row
                retry_saved += 1
                found = True; break
        if not found:
            retry_dropped += 1

    for cid, (compound, row) in assigned_ids.items():
        if compound not in final_kept:
            final_kept[compound] = row

    return list(final_kept.values()), retry_saved, retry_dropped

# ============================================================
# Main processing pipeline
# ============================================================
def process_single_mode(com_path, iden_path, score_cutoff, frag_cutoff, dccs_cutoff,
                        allowed_origins, priority_tissues, log_callback=None):
    """Process one mode (POS or NEG). Returns (filtered_rows, header_h1, header_nh, summary_dict)."""
    def log(msg):
        if log_callback: log_callback(msg)

    # Load database
    conn = open_db(); c = conn.cursor()
    ccs_hmdb = {}
    c.execute('SELECT hmdb_id, adduct, ccs, ccs_type FROM ccs WHERE hmdb_id != ""')
    for hid, ad, cv, ct in c.fetchall():
        hid = hid.upper()
        if hid not in ccs_hmdb: ccs_hmdb[hid] = []
        ccs_hmdb[hid].append((ad, cv, ct))
    hmdb_cid = {}
    c.execute('SELECT hmdb_id, pubchem_cid FROM hmdb_cid')
    for hid, pcid in c.fetchall(): hmdb_cid[hid.upper()] = pcid
    ccs_cid = {}
    c.execute('SELECT pubchem_cid, adduct, ccs, ccs_type FROM ccs WHERE pubchem_cid != ""')
    for pcid, ad, cv, ct in c.fetchall():
        if pcid not in ccs_cid: ccs_cid[pcid] = []
        ccs_cid[pcid].append((ad, cv, ct))
    ccs_name = {}
    c.execute('SELECT name, adduct, ccs, ccs_type FROM ccs WHERE name != ""')
    for nm, ad, cv, ct in c.fetchall():
        sn = simplify_name(nm)
        if sn not in ccs_name: ccs_name[sn] = []
        ccs_name[sn].append((ad, cv, ct))
    hmdb_src = {}
    c.execute('SELECT hmdb_id, origin, biospecimens, tissues FROM hmdb_source')
    for hid, o, b, t in c.fetchall(): hmdb_src[hid.upper()] = (o or '', b or '', t or '')
    conn.close()
    chemclass = load_chemclass()
    kegg = load_kegg()

    # Load COM
    t0 = time.time()
    com_ccs = {}
    com_raw = {}
    raw_headers = []
    raw_indices = []
    with open(com_path, encoding='utf-8-sig', errors='replace') as f:
        reader = csv.reader(f)
        r1 = next(reader); r2 = next(reader); r3 = next(reader)

        raw_start = None
        for i, v in enumerate(r1):
            if 'raw' in v.lower() and 'normal' not in v.lower():
                raw_start = i; break
        if raw_start is None:
            for i, v in enumerate(r2):
                if 'raw' in v.lower(): raw_start = i; break
        if raw_start is not None:
            raw_end = len(r3)
            for i, v in enumerate(r2):
                if v.strip().lower() == 'tags':
                    raw_end = i; break
            for i in range(raw_start, raw_end):
                raw_indices.append(i)
                raw_headers.append(f"Raw_{r3[i].strip()}")
        if not raw_indices:
            for i, h in enumerate(r3):
                if i >= 38 and i <= 59:
                    raw_indices.append(i)
                    raw_headers.append(f"Raw_{h.strip()}")

        for row in reader:
            if not row or not row[0].strip(): continue
            cpd = row[0].strip()
            try: com_ccs[cpd] = float(row[5].strip())
            except: pass
            com_raw[cpd] = [row[i] if i < len(row) else '' for i in raw_indices]

    has_ccs = bool(com_ccs)
    log(f"  {len(com_ccs)} CCS, {len(com_raw)} raw, {len(raw_headers)} raw cols")

    def pick_best(entries, ta):
        targets = [x.strip() for x in re.split(r'[,;]', ta) if x.strip()]
        scored = []
        for a, cv, ct in entries:
            n = norm_adduct(a); ex = 1 if ct == 'Experimental' else 0
            exact = any(norm_adduct(t) == n for t in targets)
            partial = any(n in norm_adduct(t) or norm_adduct(t) in n for t in targets)
            scored.append((ex, 2 if exact else (1 if partial else 0), cv, a, ct))
        if not scored: return None
        scored.sort(key=lambda x: (-x[0], -x[1]))
        return scored[0][2], scored[0][3], scored[0][4]

    def get_ccs(hid, name, ta):
        if not ta: return None
        if hid:
            h = hid.upper()
            entries = ccs_hmdb.get(h)
            if entries:
                x = pick_best(entries, ta)
                if x: return x
            pcid = hmdb_cid.get(h)
            if pcid:
                entries = ccs_cid.get(pcid)
                if entries:
                    x = pick_best(entries, ta)
                    if x: return x
        if name:
            entries = ccs_name.get(simplify_name(name))
            if entries:
                x = pick_best(entries, ta)
                if x: return x
        return None

    # Read IDEN
    with open(iden_path, encoding='utf-8-sig', errors='replace') as f:
        reader = csv.reader(f)
        r1 = next(reader); r2 = next(reader)
        r1_empty = sum(1 for c in r1 if not c.strip() or c.strip() in ('', '\"\"'))
        if r1_empty > len(r1) * 0.5:
            h1, h2, rows = r1, r2, list(reader)
        else:
            h1, h2, rows = None, r1, [r2] + list(reader)

    # Annotate
    total = matched = not_found = 0
    for idx, row in enumerate(rows):
        if not row or not row[0].strip(): continue
        total += 1
        compound = row[0].strip(); hid = row[1].strip()

        if has_ccs:
            ads = row[3].strip(); desc = row[13].strip(); nm = row[14].strip(); mz = row[15].strip()
            meas = com_ccs.get(compound, '')
            al = [x.strip() for x in re.split(r'[,;]', ads) if x.strip()]
            ba, _ = determine_best_adduct(nm, mz, al)
            ma = ba or ads
            res = get_ccs(hid, desc, ma)
            if res:
                tc, ta2, ct = res
                dccs = f'{abs(float(meas) - tc):.3f}' if meas else ''
                row.insert(8, str(meas) if meas else '')
                row[9] = f'{tc:.1f}'; row[10] = dccs; row.insert(11, ct)
                matched += 1
            else:
                row.insert(8, str(meas) if meas else '')
                row.insert(11, ''); not_found += 1
        else:
            row.insert(8, ''); row.insert(11, '')

        orig = bio = tis = ''
        if hid:
            s = hmdb_src.get(hid.upper())
            if s: orig, bio, tis = s
        row.insert(12, orig); row.insert(13, bio); row.insert(14, tis)
        row.insert(19, chemclass.get(hid.upper(), '') if hid else '')
        row.insert(20, kegg.get(hid.upper(), '') if hid else '')

        if (idx + 1) % 3000 == 0:
            log(f"  {idx+1}/{len(rows)} ({matched} matched, {time.time()-t0:.0f}s)")

    log(f"  Annotated: {total} rows ({time.time()-t0:.0f}s)")

    # Filter
    col_map_d = {'compound': 0, 'score': 5, 'frag': 6, 'dccs': 10,
                 'origin': 12, 'biospecimen': 13, 'tissue': 14}

    groups = {}
    for row in rows:
        if not row or not row[0].strip(): continue
        cpd = row[0].strip()
        if cpd not in groups: groups[cpd] = []
        groups[cpd].append(row)

    compound_candidates = {}
    dropped = 0
    for cpd, group_rows in groups.items():
        ranked = rank_compound_rows(group_rows, col_map_d, score_cutoff,
                                    frag_cutoff, allowed_origins, priority_tissues,
                                    dccs_cutoff)
        if ranked is None:
            dropped += 1; continue
        compound_candidates[cpd] = ranked

    log(f"  Phase A: {len(compound_candidates)}/{len(groups)} passed, {dropped} dropped")

    filtered, retry_saved, retry_dropped = hmdb_dedup(compound_candidates, col_map_d)
    log(f"  Phase B: {len(filtered)} kept [{retry_saved + retry_dropped} conflicts: {retry_saved} saved, {retry_dropped} dropped]")

    # Merge raw
    ccol = 0
    merged = 0
    for row in filtered:
        cpd = row[ccol].strip()
        ra = com_raw.get(cpd)
        if ra is not None:
            row.extend(ra)
            merged += 1
        else:
            row.extend([''] * len(raw_headers))

    missing_raw = len(filtered) - merged
    if missing_raw:
        sample = [r[ccol].strip() for r in filtered if not com_raw.get(r[ccol].strip())][:3]
        log(f"  Merged: {merged}/{len(filtered)} with raw data, {missing_raw} missing (e.g. {sample})")
    else:
        log(f"  Merged: {merged}/{len(filtered)} with raw data")

    # Build header
    nh = list(h2)
    nh.insert(8, 'Measured CCS (angstrom^2)')
    nh.insert(11, 'CCS_Type'); nh.insert(12, 'Origin')
    nh.insert(13, 'Biospecimen'); nh.insert(14, 'Tissue')
    nh.insert(19, 'Chemical Class')
    nh.insert(20, 'KEGG ID')
    nh.extend(raw_headers)

    summary = {'total': total, 'filtered': len(filtered), 'dropped': dropped,
               'retry_saved': retry_saved, 'retry_dropped': retry_dropped,
               'merged': merged, 'raw_cols': len(raw_headers)}
    return filtered, h1, nh, summary

def extract_raw_matrix(rows, header):
    """Extract Compound ID + Raw abundance matrix."""
    raw_indices = [(i, h) for i, h in enumerate(header) if h.startswith('Raw_')]
    matrix_rows = []
    id_col = 1
    for row in rows:
        cid = row[id_col] if id_col < len(row) else ''
        raw_vals = [row[i] if i < len(row) else '' for i, _ in raw_indices]
        matrix_rows.append([cid] + raw_vals)
    matrix_header = ['Compound_ID'] + [h for _, h in raw_indices]
    return matrix_header, matrix_rows
