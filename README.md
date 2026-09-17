# Virtual Embryo external-data catalog

I originally built this as a neat list of plausible outside datasets. That was the wrong standard. A useful catalog should tell you what I actually downloaded, what the files really looked like, what broke, and which sources are still only candidates.

So the repository is now **evidence-first**. Every recommended source has a committed validation manifest. Large or awkward sources stay clearly marked as unvalidated. One source, Tabula Muris, was demoted after the advertised processed-data routes failed during real validation and a more relevant developmental-heart dataset turned out to be available.

Rules snapshot: **2026-09-17**. Re-check the [official challenge rules](https://virtualembryo.ai/challenge/rules) before a final submission.

## What has actually been run

| Source | Real release tested | What happened | VEC heart-panel overlap | Spatial reality | Status |
|---|---|---|---:|---|---|
| **GSE282547 developing heart** | GSM8645799, E14.5 Visium, CD1 WT | processor passed on 456 spots × 32,285 genes | **500 / 500** | real 2-D Visium coordinates; no z invented | **recommended / validated** |
| **GSE247450 MERFISH** | GSM7890128, E9.5 WT region 0 | processor passed on 1,077 cells × 315 measured genes | **68 / 500** | metadata has x/y but no z; no `spatial_3D` invented | **validated, narrow panel** |
| **Mouse GO** | current `MOUSE-mod.gaf.gz` | processor passed; 1,111 edges for Gata4/Ctnnb1/Mab21l2 probe | n/a | n/a | **recommended / validated** |
| **sc3D E9.0** | official Figshare h5ad discovered | object is 12.19 GB, so hosted validator did not pretend to run it | not measured here | intended 3-D reconstruction | **candidate, large release** |
| **Extended Mouse Atlas** | official download index inspected | release is a 25 GB archive; not downloaded by hosted validator | not measured here | no | **candidate, large release** |
| **Tabula Muris** | old + current documented processed-data routes tried | both routes returned access failures from GitHub Actions | not measured | no | **demoted** |

The exact file URLs, SHA-256 hashes, observed schemas, commands, counts, and failures live in [`validation/`](validation/). For the two spatial expression sources, the VEC panel is read directly from the pinned public `veckit` `sample_heart_9.5.h5ad` rather than from a hand-copied gene list.

That table is the main product of this repo. The processors are useful because they make those claims reproducible.

## The sources I would actually start with

### 1. GSE282547 developing mouse-heart atlas

Source: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282547>

This became the most interesting expression/spatial source after actually testing it. The associated 2026 study profiled developing mouse hearts with Visium and Slide-seq. For VEC, the later allowed stages are much more relevant than an adult atlas.

The validated sample is **GSM8645799, E14.5 Visium, CD1 wild type**. The real release contained:

- 456 spots
- 32,285 genes
- all **500/500** genes in the pinned public VEC heart panel
- 2-D pixel coordinates from the Visium tissue-position file
- no 3-D coordinate field

The adapter deliberately stores those coordinates as `obsm['spatial_2D']`. It does not manufacture a z-axis just because VEC Task 2 submissions are 3-D.

```bash
python processors/process_gse282547_visium.py raw/GSM8645799 \
  --stage 14.5 \
  --task t2-heart \
  --gene-list vec_heart_500.txt \
  --out processed/GSM8645799.E14.5.h5ad
```

Why I like it: it is heart-specific, developmental, later than the challenge's protected extrapolation window, and genome-wide enough to contain the whole public heart panel.

What it does **not** give you: a ready-made 3-D embryo. It is a 2-D section, so using it for T2 geometry requires an actual modeling idea rather than a schema conversion trick.

Evidence: [`validation/gse282547_e14_5_visium.json`](validation/gse282547_e14_5_visium.json)

### 2. GSE247450 embryonic endothelial MERFISH

Source: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247450>

This still earns a place because it is embryonic MERFISH and therefore closer to the assay style of the spatial challenge data. Real validation changed how I would describe its usefulness, though.

For GSM7890128 E9.5 region 0, the released matrix has **315 genes**, and only **68/500** overlap the pinned public VEC heart panel. The metadata has 22 observed columns including `center_x` and `center_y`, but no genuine z coordinate.

So I would treat this as a narrow same-modality / endothelial prior, not as a drop-in 500-gene training set.

```bash
python processors/process_gse247450.py raw/GSE247450 \
  --task t2-heart \
  --gene-list vec_heart_500.txt \
  --out-dir processed/gse247450
```

The adapter pairs the real matrix/metadata files by sample key, checks cell IDs to infer orientation, and leaves `spatial_3D` absent when the source does not provide z.

Evidence: [`validation/gse247450_e9_5_region0.json`](validation/gse247450_e9_5_region0.json)

### 3. Mouse Gene Ontology annotations

Current mouse GAF: <https://current.geneontology.org/annotations/gaf/MOUSE-mod.gaf.gz>

This is useful for a different reason. Task 3 provides extremely little supervision over perturbation identity. GO gives gene-function structure without pretending to be another embryo measurement.

The real validation run downloaded and hashed the current GAF, then tested `Gata4`, `Ctnnb1`, and `Mab21l2`. It produced **1,111 gene-GO edges** covering all three probe genes.

```bash
python processors/prepare_mouse_go.py MOUSE-mod.gaf.gz \
  --gene-list vec_500_genes.txt \
  --exclude-iea \
  --out processed/mouse_go_vec500.csv
```

Keep the evidence codes. Whether to exclude IEA annotations is a modeling choice.

Evidence: [`validation/mouse_go.json`](validation/mouse_go.json)

## Good candidates I have not fully run

### Extended Mouse Atlas

Project: <https://marionilab.github.io/ExtendedMouseAtlas/>

The scientific fit is strong: 430,339 cells over 13 time points from E6.5 through E9.5. The practical problem is that the official downloadable package is a **25 GB archive**. The processor has synthetic end-to-end coverage, but I have not run it against that full archive in the hosted validator.

```bash
python processors/process_extended_mouse_atlas.py embryo_complete.h5ad \
  --task t1 \
  --out processed/extended_mouse_atlas.t1.h5ad
```

I would use it before Tabula Muris if I were willing to handle the download. I would not describe this adapter as real-release validated yet.

Evidence: [`validation/extended_mouse_atlas.json`](validation/extended_mouse_atlas.json)

### sc3D / GSE197353

GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197353>

The convenient E9.0 reconstructed h5ad exposed through Figshare is **12,187,200,243 bytes**. That is large enough that the hosted validator records the exact release metadata instead of downloading it and pretending a compile test proves the adapter.

The source remains scientifically attractive for 3-D spatial priors. If I used it in a serious model, the next step would be running the processor on the real object on a machine with enough disk/RAM and committing the resulting manifest.

Evidence: [`validation/sc3d_e9_0.json`](validation/sc3d_e9_0.json)

## What I dropped down the list

### Tabula Muris

I initially included Tabula Muris because it gave the catalog a broad adult tissue prior. That was taxonomy-first curation.

In practice it is developmentally far from VEC, and both processed-data access routes I tried failed from the hosted validator:

- legacy `czbiohub-tabula-muris` route: HTTP/S3 403
- current Open Data Registry `czb-tabula-muris` object route: `HeadObject 403 Forbidden`

I stopped there because GSE282547 is substantially more task-relevant. The old processor remains for anyone who already has the files, but `audit_catalog.py` hides this entry by default.

Evidence: [`validation/tabula_muris.json`](validation/tabula_muris.json)

## Safety behavior

The preprocessing helpers are intentionally conservative.

- Ambiguous stage strings such as ranges are rejected instead of silently taking the first number.
- Somite/Theiler staging needs a source-specific parser.
- T1/T2 stage windows use explicit open/closed interval logic with regression tests around the boundaries.
- Task 3 cannot be declared safe from stage alone. T3 processor runs emit `T3_GENOTYPE_NOT_AUDITED` and record that status in provenance.
- A source with x/y coordinates does not automatically become 3-D. Adapters only create `spatial_3D` when a real 3-D field has been validated.

The rules change independently of this repo. Treat the challenge site as authoritative.

## Two kinds of testing

### Synthetic CI

Every processor has a small end-to-end fixture. CI installs the real runtime dependencies and actually reads/writes AnnData or source-format files. These tests catch parser, filtering, schema, coordinate, and provenance regressions without downloading multi-GB releases.

```bash
pytest -q tests
```

### Real-source validation

[`scripts/validate_real_sources.py`](scripts/validate_real_sources.py) downloads specific public releases where practical, executes the processor, and commits observed evidence. [`scripts/add_vec_panel_overlap.py`](scripts/add_vec_panel_overlap.py) separately compares real source genes against the pinned public VEC heart panel.

A validation failure stays visible. The workflow commits whatever evidence it obtained and then fails, rather than turning a broken source into a missing log.

Typical manifest fields include:

```json
{
  "tested_on": "2026-09-17",
  "status": "real_release_processor_passed",
  "source_files": [{"filename": "...", "sha256": "...", "bytes": 123}],
  "output_cells": 1077,
  "observed_metadata_columns": ["..."],
  "vec_heart_panel": {"panel_genes": 500, "genes_found": 68}
}
```

## Source terms

I resolved the easy cases instead of repeating “check the licence” everywhere. See [`SOURCE_TERMS.md`](SOURCE_TERMS.md) for the exact status and links.

The short version:

- NCBI itself places no restrictions on use/distribution of GEO molecular data, while warning that submitters may retain third-party rights.
- Gene Ontology data products are CC BY 4.0 and should be attributed with the release/date where possible.
- I did not find a dataset-specific redistribution licence on the Extended Mouse Atlas data page, so this repo links upstream rather than vendoring the atlas.
- The GSE282547 paper's CC BY-NC-ND licence is an **article** licence. I do not treat it as a blanket licence for every separately deposited GEO file.

## Machine-readable catalog

[`catalog.json`](catalog.json) records the evidence tier, validation manifest, task notes, processor, and source terms for each entry.

```bash
python audit_catalog.py --task t2-heart
```

Demoted entries are hidden unless you explicitly ask for them:

```bash
python audit_catalog.py --task t2-heart --include-demoted
```

## Scope

This is deliberately not an exhaustive atlas directory. I would rather have three sources whose quirks are documented than thirty links that merely look relevant.

An entry being **validated** means the advertised processor was exercised on the named real release and the evidence was committed. It does not mean the source will improve a model, that every file from the study is safe to use, or that the challenge organizers have pre-approved it.
