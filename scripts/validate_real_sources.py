#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
VALIDATION = ROOT / "validation"
VALIDATION.mkdir(exist_ok=True)
WORK = Path(tempfile.mkdtemp(prefix="vec-external-real-"))
TESTED_ON = datetime.now(timezone.utc).date().isoformat()


def fetch(url: str, dest: Path) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    h = hashlib.sha256()
    total = 0
    req = urllib.request.Request(url, headers={"User-Agent": "vec-external-catalog-validation/1.0"})
    with urllib.request.urlopen(req, timeout=120) as r, dest.open("wb") as f:
        while True:
            chunk = r.read(8 * 1024 * 1024)
            if not chunk:
                break
            f.write(chunk)
            h.update(chunk)
            total += len(chunk)
    return {"url": url, "sha256": h.hexdigest(), "bytes": total, "filename": dest.name}


def head_size(url: str) -> int | None:
    try:
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "vec-external-catalog-validation/1.0"})
        with urllib.request.urlopen(req, timeout=60) as r:
            value = r.headers.get("Content-Length")
            return int(value) if value else None
    except Exception:
        return None


def run(*args: object) -> subprocess.CompletedProcess:
    cmd = [sys.executable, *map(str, args)]
    return subprocess.run(cmd, cwd=ROOT, check=True, text=True, capture_output=True)


def write_manifest(name: str, payload: dict) -> None:
    payload = {"tested_on": TESTED_ON, **payload}
    (VALIDATION / f"{name}.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def validate_gse247450() -> None:
    raw = WORK / "gse247450"
    urls = [
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7890nnn/GSM7890128/suppl/GSM7890128_E9.5_region_0_cell_by_gene.csv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7890nnn/GSM7890128/suppl/GSM7890128_E9.5_region_0_cell_metadata.csv.gz",
    ]
    files = [fetch(url, raw / url.rsplit("/", 1)[-1]) for url in urls]
    matrix_path = raw / files[0]["filename"]
    meta_path = raw / files[1]["filename"]
    matrix = pd.read_csv(matrix_path, index_col=0)
    meta = pd.read_csv(meta_path, index_col=0)

    out_dir = WORK / "gse247450_out"
    proc = run(ROOT / "processors/process_gse247450.py", raw, "--task", "t2-heart", "--out-dir", out_dir)
    out = out_dir / "GSM7890128_E9.5_region_0.t2-heart.h5ad"
    result = ad.read_h5ad(out)
    provenance = json.loads(out.with_suffix(".provenance.json").read_text())

    write_manifest("gse247450_e9_5_region0", {
        "source": "GSE247450 / GSM7890128",
        "source_files": files,
        "command": "python processors/process_gse247450.py <raw_dir> --task t2-heart --out-dir <out_dir>",
        "processor_exit": proc.returncode,
        "processor_stdout": proc.stdout.strip(),
        "raw_matrix_shape_as_stored": list(matrix.shape),
        "raw_metadata_shape": list(meta.shape),
        "observed_metadata_columns": list(map(str, meta.columns)),
        "output_cells": int(result.n_obs),
        "output_genes": int(result.n_vars),
        "spatial_3d_present": "spatial_3D" in result.obsm,
        "provenance": provenance,
    })


def validate_mouse_go() -> None:
    url = "https://current.geneontology.org/annotations/gaf/MOUSE-mod.gaf.gz"
    src = WORK / "MOUSE-mod.gaf.gz"
    source = fetch(url, src)
    genes = WORK / "go_probe_genes.txt"
    genes.write_text("Gata4\nCtnnb1\nMab21l2\n", encoding="utf-8")
    out = WORK / "mouse_go_probe.csv"
    proc = run(ROOT / "processors/prepare_mouse_go.py", src, "--gene-list", genes, "--out", out)
    table = pd.read_csv(out)
    write_manifest("mouse_go", {
        "source": "Mouse Gene Ontology GAF",
        "source_files": [source],
        "command": "python processors/prepare_mouse_go.py MOUSE-mod.gaf.gz --gene-list <3-gene probe> --out <csv>",
        "probe_genes": ["Gata4", "Ctnnb1", "Mab21l2"],
        "processor_exit": proc.returncode,
        "processor_stdout": proc.stdout.strip(),
        "output_rows": int(len(table)),
        "genes_found": sorted(table["gene"].astype(str).unique().tolist()) if len(table) else [],
        "evidence_codes_observed": sorted(table["evidence"].astype(str).unique().tolist()) if len(table) else [],
    })


def validate_tabula_muris() -> None:
    matrix_url = "https://s3.amazonaws.com/czbiohub-tabula-muris/TM_droplet_mat.h5ad"
    meta_url = "https://s3.amazonaws.com/czbiohub-tabula-muris/TM_droplet_metadata.csv"
    advertised_size = head_size(matrix_url)
    # This is deliberately a real-run gate, not a fake success. If the release becomes too
    # large for the hosted validation runner, record that fact and stop advertising it as validated.
    if advertised_size is not None and advertised_size > 2_500_000_000:
        write_manifest("tabula_muris", {
            "source": "Tabula Muris droplet release",
            "status": "not_run_release_too_large_for_hosted_validator",
            "matrix_url": matrix_url,
            "content_length": advertised_size,
        })
        return

    matrix = WORK / "TM_droplet_mat.h5ad"
    metadata = WORK / "TM_droplet_metadata.csv"
    sources = [fetch(matrix_url, matrix), fetch(meta_url, metadata)]
    meta = pd.read_csv(metadata, index_col=0)
    tissue_col = "tissue" if "tissue" in meta.columns else None
    if tissue_col is None:
        raise RuntimeError(f"Expected tissue column in Tabula Muris metadata; saw {list(meta.columns)}")
    tissue_values = sorted({str(x) for x in meta[tissue_col].dropna().unique()})
    heart_values = [x for x in tissue_values if "heart" in x.lower()]
    if not heart_values:
        raise RuntimeError(f"No heart-like tissue label in real metadata: {tissue_values}")

    out = WORK / "tabula_muris_heart.h5ad"
    cmd = [ROOT / "processors/process_tabula_muris.py", matrix, "--metadata", metadata]
    for tissue in heart_values:
        cmd += ["--tissue", tissue]
    cmd += ["--out", out]
    proc = run(*cmd)
    inp = ad.read_h5ad(matrix, backed="r")
    result = ad.read_h5ad(out)
    provenance = json.loads(out.with_suffix(".provenance.json").read_text())
    write_manifest("tabula_muris", {
        "source": "Tabula Muris droplet release",
        "source_files": sources,
        "command": "python processors/process_tabula_muris.py TM_droplet_mat.h5ad --metadata TM_droplet_metadata.csv --tissue <observed heart label(s)> --out <h5ad>",
        "processor_exit": proc.returncode,
        "processor_stdout": proc.stdout.strip(),
        "input_cells": int(inp.n_obs),
        "input_genes": int(inp.n_vars),
        "metadata_rows": int(len(meta)),
        "observed_metadata_columns": list(map(str, meta.columns)),
        "heart_labels_used": heart_values,
        "output_cells": int(result.n_obs),
        "output_genes": int(result.n_vars),
        "provenance": provenance,
    })


def probe_sc3d_release() -> None:
    api = "https://api.figshare.com/v2/articles/21695879/versions/1/files"
    req = urllib.request.Request(api, headers={"User-Agent": "vec-external-catalog-validation/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        files = json.load(r)
    h5ads = [f for f in files if str(f.get("name", "")).lower().endswith(".h5ad")]
    if not h5ads:
        write_manifest("sc3d_e9_0", {"source": "sc3D E9.0 Figshare", "status": "no_h5ad_found_via_figshare_api", "api": api, "files": files})
        return
    item = h5ads[0]
    size = int(item.get("size") or 0)
    # Try the real release when it is reasonable for the hosted runner. Otherwise commit a
    # truthful size/status manifest rather than implying that the adapter was exercised.
    if size > 2_500_000_000:
        write_manifest("sc3d_e9_0", {
            "source": "sc3D E9.0 Figshare",
            "status": "not_run_release_too_large_for_hosted_validator",
            "figshare_file": item,
        })
        return
    src = WORK / str(item["name"])
    source = fetch(str(item["download_url"]), src)
    out_dir = WORK / "sc3d_out"
    proc = run(ROOT / "processors/process_sc3d.py", src, "--task", "t2-heart", "--stage", "9.0", "--out-dir", out_dir)
    outputs = list(out_dir.glob("*.h5ad"))
    if len(outputs) != 1:
        raise RuntimeError(f"Expected one sc3D output, got {outputs}")
    result = ad.read_h5ad(outputs[0])
    provenance = json.loads(outputs[0].with_suffix(".provenance.json").read_text())
    write_manifest("sc3d_e9_0", {
        "source": "sc3D E9.0 Figshare",
        "figshare_file": item,
        "source_files": [source],
        "command": "python processors/process_sc3d.py <E9.0 h5ad> --task t2-heart --stage 9.0 --out-dir <out_dir>",
        "processor_exit": proc.returncode,
        "processor_stdout": proc.stdout.strip(),
        "output_cells": int(result.n_obs),
        "output_genes": int(result.n_vars),
        "spatial_3d_present": "spatial_3D" in result.obsm,
        "provenance": provenance,
    })


def record_extended_mouse_atlas_constraint() -> None:
    # The official release page distributes the package as a 25 GB tarball. Do not call the
    # adapter "real-run validated" until we have actually processed that release elsewhere.
    write_manifest("extended_mouse_atlas", {
        "source": "Extended Mouse Atlas",
        "status": "not_real_run_validated",
        "reason": "official downloadable package is a 25 GB tarball; not downloaded by the hosted validator",
        "release_index": "https://bioinformatics.stemcells.cam.ac.uk/rlh60/Supplemental/ExtendedMouseAtlas/",
        "advertised_archive": "ExtendedMouseAtlas.tar.gz",
        "advertised_size": "25G",
    })


def main() -> None:
    print("workdir", WORK)
    validate_gse247450()
    validate_mouse_go()
    validate_tabula_muris()
    probe_sc3d_release()
    record_extended_mouse_atlas_constraint()
    print("wrote manifests:")
    for p in sorted(VALIDATION.glob("*.json")):
        print(" -", p.relative_to(ROOT))


if __name__ == "__main__":
    main()
