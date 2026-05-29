"""Metabolomics single-omics analysis helpers for the web app."""
import csv
import itertools
import math
import os
import re
from collections import Counter

import numpy as np
import pandas as pd


ID_COLUMN_CANDIDATES = ("Compound_ID", "Compound ID", "HMDB_ID", "HMDB ID", "ID")


def parse_sample_name(name):
    """Parse sample metadata from names such as Raw_20250528-HXY-NEG-C1-1."""
    clean = str(name).strip()
    meta = {
        "sample_id": clean,
        "date": "",
        "project": "",
        "mode": "",
        "group": "UNPARSED",
        "replicate": "",
        "is_qc": False,
    }
    if clean.startswith("Raw_"):
        parts = clean[4:].split("-")
        if len(parts) >= 5:
            meta["date"] = parts[0]
            meta["project"] = parts[1]
            meta["mode"] = parts[2]
            meta["group"] = "-".join(parts[3:-1]) or "UNPARSED"
            meta["replicate"] = parts[-1]
    else:
        match = re.search(r"(.+?)[._-](\d+)$", clean)
        if match:
            meta["group"] = match.group(1)
            meta["replicate"] = match.group(2)

    meta["is_qc"] = meta["group"].upper().startswith("QC")
    return meta


def _find_id_column(columns):
    for candidate in ID_COLUMN_CANDIDATES:
        if candidate in columns:
            return candidate
    return columns[0]


def _detect_file_mode(filename):
    base = os.path.basename(filename or "").upper()
    if re.search(r"(^|[_\-.])POS([_\-.]|$)", base):
        return "POS"
    if re.search(r"(^|[_\-.])NEG([_\-.]|$)", base):
        return "NEG"
    return ""


def read_matrix(path):
    df = pd.read_csv(path, encoding="utf-8-sig")
    if df.empty:
        raise ValueError("The matrix is empty.")
    id_col = _find_id_column(list(df.columns))
    sample_cols = [c for c in df.columns if c != id_col]
    if not sample_cols:
        raise ValueError("No sample columns were found.")

    values = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    ids = df[id_col].astype(str).fillna("")
    samples = [parse_sample_name(c) for c in sample_cols]
    return df, ids, values, id_col, sample_cols, samples


def _safe_float(value, digits=4):
    if value is None:
        return None
    try:
        if not np.isfinite(value):
            return None
        return round(float(value), digits)
    except Exception:
        return None


def _bh_fdr(p_values):
    p = np.asarray(p_values, dtype=float)
    q = np.ones_like(p)
    finite = np.isfinite(p)
    if not finite.any():
        return q
    idx = np.where(finite)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order]
    m = len(ranked)
    adjusted = np.empty(m, dtype=float)
    running = 1.0
    for i in range(m - 1, -1, -1):
        running = min(running, ranked[i] * m / (i + 1))
        adjusted[i] = running
    q[order] = np.clip(adjusted, 0, 1)
    q[~finite] = 1
    return q


def _welch_stat(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return 0.0
    va = float(np.var(a, ddof=1))
    vb = float(np.var(b, ddof=1))
    denom = math.sqrt(va / len(a) + vb / len(b))
    diff = abs(float(np.mean(a) - np.mean(b)))
    if denom == 0:
        return 0.0 if diff == 0 else float("inf")
    return diff / denom


def _permutation_pvalue(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    pooled = np.concatenate([a[np.isfinite(a)], b[np.isfinite(b)]])
    n_a = len(a[np.isfinite(a)])
    n_b = len(b[np.isfinite(b)])
    if n_a < 2 or n_b < 2:
        return 1.0
    total = math.comb(n_a + n_b, n_a)
    obs = _welch_stat(pooled[:n_a], pooled[n_a:])

    if total <= 5000:
        count = 0
        for idx in itertools.combinations(range(n_a + n_b), n_a):
            mask = np.zeros(n_a + n_b, dtype=bool)
            mask[list(idx)] = True
            stat = _welch_stat(pooled[mask], pooled[~mask])
            if stat >= obs - 1e-12:
                count += 1
        return count / total

    rng = np.random.default_rng(42)
    count = 0
    rounds = 5000
    for _ in range(rounds):
        perm = rng.permutation(n_a + n_b)
        stat = _welch_stat(pooled[perm[:n_a]], pooled[perm[n_a:]])
        if stat >= obs - 1e-12:
            count += 1
    return (count + 1) / (rounds + 1)


def _normalize_for_pca(values):
    x = np.log2(values.fillna(0).clip(lower=0).to_numpy(dtype=float) + 1.0)
    keep = np.nanvar(x, axis=1) > 0
    x = x[keep, :]
    if x.shape[0] == 0:
        return None
    x = x.T
    x = x - np.nanmean(x, axis=0)
    sd = np.nanstd(x, axis=0, ddof=0)
    sd[sd == 0] = 1
    x = np.nan_to_num(x / sd)
    return x


def pca_scores(values, sample_cols, samples, use_biological=True):
    selected = [
        i for i, sample in enumerate(samples)
        if (not use_biological or not sample["is_qc"])
    ]
    if len(selected) < 2:
        return {"points": [], "explained": [0, 0]}

    selected_cols = [sample_cols[i] for i in selected]
    x = _normalize_for_pca(values[selected_cols])
    if x is None or x.shape[0] < 2:
        return {"points": [], "explained": [0, 0]}

    _, s, vt = np.linalg.svd(x, full_matrices=False)
    scores = np.dot(x, vt.T[:, :2])
    denom = np.sum(s ** 2)
    explained = (s[:2] ** 2 / denom * 100) if denom else np.array([0, 0])
    if len(explained) == 1:
        explained = np.array([explained[0], 0])

    points = []
    for row_idx, sample_idx in enumerate(selected):
        meta = samples[sample_idx]
        points.append({
            "sample": sample_cols[sample_idx],
            "group": meta["group"],
            "replicate": meta["replicate"],
            "pc1": _safe_float(scores[row_idx, 0]),
            "pc2": _safe_float(scores[row_idx, 1] if scores.shape[1] > 1 else 0),
        })
    return {
        "points": points,
        "explained": [_safe_float(explained[0], 2), _safe_float(explained[1], 2)],
    }


def qc_rsd(values, samples):
    qc_cols = [s["sample_id"] for s in samples if s["group"].upper() == "QC"]
    if len(qc_cols) < 2:
        return pd.Series([np.nan] * values.shape[0], index=values.index)
    qc = values[qc_cols].replace(0, np.nan)
    mean = qc.mean(axis=1, skipna=True)
    sd = qc.std(axis=1, skipna=True)
    return sd / mean * 100


def summarize(path, ids, values, sample_cols, samples):
    groups = Counter(s["group"] for s in samples)
    modes = sorted({s["mode"] for s in samples if s["mode"]})
    file_mode = _detect_file_mode(path)
    warnings = []
    if file_mode and modes and any(mode.upper() != file_mode for mode in modes):
        warnings.append(
            f"File name suggests {file_mode}, but sample names contain mode(s): {', '.join(modes)}."
        )
    if any(s["group"] == "UNPARSED" for s in samples):
        warnings.append("Some sample names could not be parsed into group and replicate.")

    numeric_total = values.size
    missing = int(values.isna().sum().sum())
    zero = int((values.fillna(0) == 0).sum().sum())
    biological = [s for s in samples if not s["is_qc"]]
    qc = [s for s in samples if s["is_qc"]]

    sample_table = []
    for s in samples:
        sample_table.append({
            "sample_id": s["sample_id"],
            "date": s["date"],
            "project": s["project"],
            "mode": s["mode"],
            "group": s["group"],
            "replicate": s["replicate"],
            "is_qc": s["is_qc"],
        })

    return {
        "feature_count": int(len(ids)),
        "sample_count": int(len(sample_cols)),
        "biological_sample_count": int(len(biological)),
        "qc_sample_count": int(len(qc)),
        "missing_percent": _safe_float(missing / numeric_total * 100, 2),
        "zero_percent": _safe_float(zero / numeric_total * 100, 2),
        "groups": [{"group": k, "count": groups[k]} for k in sorted(groups)],
        "modes": modes,
        "file_mode": file_mode,
        "warnings": warnings,
        "samples": sample_table,
    }


def differential(values, ids, samples, group_a, group_b):
    group_to_cols = {}
    for s in samples:
        if not s["is_qc"]:
            group_to_cols.setdefault(s["group"], []).append(s["sample_id"])
    if group_a not in group_to_cols or group_b not in group_to_cols:
        raise ValueError("Selected comparison groups were not found in biological samples.")

    cols_a = group_to_cols[group_a]
    cols_b = group_to_cols[group_b]
    a_raw = values[cols_a].fillna(0).clip(lower=0)
    b_raw = values[cols_b].fillna(0).clip(lower=0)
    a_log = np.log2(a_raw.to_numpy(dtype=float) + 1.0)
    b_log = np.log2(b_raw.to_numpy(dtype=float) + 1.0)

    mean_a_raw = a_raw.mean(axis=1).to_numpy(dtype=float)
    mean_b_raw = b_raw.mean(axis=1).to_numpy(dtype=float)
    log2_fc = np.log2((mean_b_raw + 1.0) / (mean_a_raw + 1.0))
    p_values = np.array([_permutation_pvalue(a_log[i, :], b_log[i, :]) for i in range(values.shape[0])])
    q_values = _bh_fdr(p_values)

    result = pd.DataFrame({
        "Compound_ID": ids.to_numpy(),
        f"mean_{group_a}": mean_a_raw,
        f"mean_{group_b}": mean_b_raw,
        "log2FC": log2_fc,
        "p_value": p_values,
        "FDR": q_values,
    })
    result["direction"] = np.where(result["log2FC"] > 0, "up", np.where(result["log2FC"] < 0, "down", "flat"))
    result["abs_log2FC"] = result["log2FC"].abs()
    result = result.sort_values(["FDR", "p_value", "abs_log2FC"], ascending=[True, True, False])
    return result


def heatmap_payload(values, ids, samples, diff_result, max_features=30):
    biological_cols = [s["sample_id"] for s in samples if not s["is_qc"]]
    if not biological_cols:
        return {"features": [], "samples": [], "groups": [], "matrix": []}
    top_ids = diff_result.head(max_features)["Compound_ID"].astype(str).tolist()
    id_to_pos = {str(cid): i for i, cid in enumerate(ids)}
    selected_pos = [id_to_pos[cid] for cid in top_ids if cid in id_to_pos]
    if not selected_pos:
        return {"features": [], "samples": biological_cols, "groups": [], "matrix": []}
    x = np.log2(values.iloc[selected_pos][biological_cols].fillna(0).clip(lower=0).to_numpy(dtype=float) + 1.0)
    row_mean = np.mean(x, axis=1, keepdims=True)
    row_sd = np.std(x, axis=1, keepdims=True)
    row_sd[row_sd == 0] = 1
    z = np.clip((x - row_mean) / row_sd, -3, 3)
    group_map = {s["sample_id"]: s["group"] for s in samples}
    return {
        "features": [str(ids.iloc[i]) for i in selected_pos],
        "samples": biological_cols,
        "groups": [group_map[c] for c in biological_cols],
        "matrix": [[_safe_float(v, 3) for v in row] for row in z],
    }


def volcano_payload(diff_result):
    rows = []
    for _, row in diff_result.iterrows():
        p = max(float(row["p_value"]), 1e-300)
        rows.append({
            "id": str(row["Compound_ID"]),
            "log2FC": _safe_float(row["log2FC"], 4),
            "negLog10P": _safe_float(-math.log10(p), 4),
            "FDR": _safe_float(row["FDR"], 4),
            "significant": bool(row["FDR"] <= 0.05 and abs(row["log2FC"]) >= 1),
        })
    return rows


def _plot_pca(pca, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    points = pca.get("points", [])
    if not points:
        return
    groups = sorted({p["group"] for p in points})
    cmap = plt.get_cmap("tab10")
    colors = {g: cmap(i % 10) for i, g in enumerate(groups)}
    fig, ax = plt.subplots(figsize=(6.2, 4.6), dpi=180)
    for g in groups:
        pts = [p for p in points if p["group"] == g]
        ax.scatter([p["pc1"] for p in pts], [p["pc2"] for p in pts], s=42, label=g, color=colors[g])
        for p in pts:
            ax.annotate(p["replicate"], (p["pc1"], p["pc2"]), fontsize=7, xytext=(3, 3), textcoords="offset points")
    explained = pca.get("explained", [0, 0])
    ax.set_xlabel(f"PC1 ({explained[0]}%)")
    ax.set_ylabel(f"PC2 ({explained[1]}%)")
    ax.axhline(0, color="#d0d7de", linewidth=0.8)
    ax.axvline(0, color="#d0d7de", linewidth=0.8)
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.set_title("PCA")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _plot_volcano(volcano, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not volcano:
        return
    x = [r["log2FC"] for r in volcano]
    y = [r["negLog10P"] for r in volcano]
    colors = []
    for r in volcano:
        if r["significant"] and r["log2FC"] > 0:
            colors.append("#dc2626")
        elif r["significant"] and r["log2FC"] < 0:
            colors.append("#2563eb")
        else:
            colors.append("#94a3b8")
    fig, ax = plt.subplots(figsize=(6.2, 4.6), dpi=180)
    ax.scatter(x, y, s=14, c=colors, alpha=0.78, edgecolors="none")
    ax.axvline(-1, color="#d0d7de", linestyle="--", linewidth=0.8)
    ax.axvline(1, color="#d0d7de", linestyle="--", linewidth=0.8)
    ax.axhline(-math.log10(0.05), color="#d0d7de", linestyle="--", linewidth=0.8)
    ax.set_xlabel("log2FC")
    ax.set_ylabel("-log10(P)")
    ax.set_title("Volcano")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _plot_heatmap(heatmap, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matrix = heatmap.get("matrix", [])
    if not matrix:
        return
    data = np.asarray(matrix, dtype=float)
    fig_h = max(4.5, min(10, data.shape[0] * 0.22 + 1.6))
    fig, ax = plt.subplots(figsize=(8.6, fig_h), dpi=180)
    im = ax.imshow(data, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
    def short_label(sample):
        parts = str(sample).split("-")
        if len(parts) >= 2:
            return parts[-2] + "-" + parts[-1]
        return str(sample)

    ax.set_xticks(range(len(heatmap["samples"])))
    ax.set_xticklabels([short_label(s) for s in heatmap["samples"]], rotation=60, ha="right", fontsize=7)
    ax.set_yticks(range(len(heatmap["features"])))
    ax.set_yticklabels(heatmap["features"], fontsize=7)
    ax.set_title("Top Differential Feature Heatmap")
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label="Z-score")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _write_report(output_dir, payload):
    report_path = os.path.join(output_dir, "analysis_report.txt")
    summary = payload["summary"]
    comparison = payload["comparison"]
    lines = [
        "OmicsBridge Metabolomics Analysis Report",
        "=" * 42,
        f"Features: {summary['feature_count']}",
        f"Samples: {summary['sample_count']}",
        f"Biological samples: {summary['biological_sample_count']}",
        f"QC samples: {summary['qc_sample_count']}",
        f"Zero percent: {summary['zero_percent']}%",
        f"Missing percent: {summary['missing_percent']}%",
        f"Comparison: {comparison['group_a']} vs {comparison['group_b']}",
        f"Significant features: {payload['significant_count']} (FDR <= 0.05 and |log2FC| >= 1)",
        "",
        "Groups:",
    ]
    lines.extend([f"  {g['group']}: {g['count']}" for g in summary["groups"]])
    if summary["warnings"]:
        lines.extend(["", "Warnings:"])
        lines.extend([f"  {w}" for w in summary["warnings"]])
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return report_path


def run_analysis(path, output_dir, group_a=None, group_b=None, write_figures=True):
    _, ids, values, id_col, sample_cols, samples = read_matrix(path)
    summary = summarize(path, ids, values, sample_cols, samples)
    biological_groups = [g["group"] for g in summary["groups"] if not g["group"].upper().startswith("QC")]
    if not group_a or not group_b:
        if len(biological_groups) < 2:
            raise ValueError("At least two biological groups are required for differential analysis.")
        group_a, group_b = biological_groups[0], biological_groups[1]

    rsd = qc_rsd(values, samples)
    diff = differential(values, ids, samples, group_a, group_b)
    diff.insert(1, "QC_RSD_percent", rsd.loc[diff.index].to_numpy())

    os.makedirs(output_dir, exist_ok=True)
    diff_path = os.path.join(output_dir, f"differential_{group_a}_vs_{group_b}.csv")
    sample_path = os.path.join(output_dir, "sample_design.csv")
    qc_path = os.path.join(output_dir, "qc_summary.csv")

    diff.to_csv(diff_path, index=False, encoding="utf-8-sig")
    with open(sample_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=["sample_id", "date", "project", "mode", "group", "replicate", "is_qc"])
        writer.writeheader()
        writer.writerows(summary["samples"])

    qc_df = pd.DataFrame({
        "Compound_ID": ids,
        "missing_percent": values.isna().mean(axis=1) * 100,
        "zero_percent": (values.fillna(0) == 0).mean(axis=1) * 100,
        "QC_RSD_percent": rsd,
    })
    qc_df.to_csv(qc_path, index=False, encoding="utf-8-sig")

    sig = diff[(diff["FDR"] <= 0.05) & (diff["abs_log2FC"] >= 1)]
    payload = {
        "id_column": id_col,
        "summary": summary,
        "comparison": {"group_a": group_a, "group_b": group_b},
        "available_groups": biological_groups,
        "qc": {
            "features_with_qc_rsd_le_30": int((rsd <= 30).sum()) if rsd.notna().any() else None,
            "features_with_qc_rsd_le_50": int((rsd <= 50).sum()) if rsd.notna().any() else None,
        },
        "pca": pca_scores(values, sample_cols, samples, use_biological=True),
        "volcano": volcano_payload(diff),
        "heatmap": heatmap_payload(values, ids, samples, diff),
        "top_table": diff.head(25).replace({np.nan: None}).to_dict(orient="records"),
        "significant_count": int(len(sig)),
        "files": [os.path.basename(diff_path), os.path.basename(sample_path), os.path.basename(qc_path)],
    }
    report_path = _write_report(output_dir, payload)
    payload["files"].append(os.path.basename(report_path))
    if write_figures:
        figure_paths = [
            ("pca.png", _plot_pca, payload["pca"]),
            ("volcano.png", _plot_volcano, payload["volcano"]),
            ("heatmap.png", _plot_heatmap, payload["heatmap"]),
        ]
        for filename, plotter, data in figure_paths:
            figure_path = os.path.join(output_dir, filename)
            plotter(data, figure_path)
            if os.path.exists(figure_path):
                payload["files"].append(filename)
    return payload
