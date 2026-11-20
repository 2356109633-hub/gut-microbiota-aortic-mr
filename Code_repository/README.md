# Gut microbiota and aortic aneurysm/dissection: a two-sample MR with independent replication

This repository contains the analysis code and summary result tables for the manuscript:

> **Gut microbiota and aortic aneurysm/dissection: a two-sample Mendelian randomization study with independent replication**

**Code DOI:** `[Zenodo DOI to be inserted after deposit]`

## Study overview

| Item | Description |
|---|---|
| Exposure | MiBioGen gut microbiome GWAS (n = 18,340; 140 bacterial taxa; 1,748 instrument variants; 1,363 unique SNPs) |
| Discovery | FinnGen release 13 (DF13), four aortic outcomes |
| Replication | AAAgen abdominal aortic aneurysm GWAS (39,221 cases; 1,086,107 controls) |
| Primary method | Inverse-variance weighted (IVW) MR |
| Sensitivity | Weighted median, MR-Egger, Cochran's Q, MR-Egger intercept, MR-PRESSO, leave-one-out, Steiger filtering |
| Software | R 4.3.1; TwoSampleMR 0.4.26; MRPRESSO 1.0 |

### Discovery outcomes (FinnGen R13)

| Outcome | Cases | Controls | Endpoint |
|---|---:|---:|---|
| Abdominal aortic aneurysm (AAA) | 4,634 | 460,490 | I9_ABAORTANEUR |
| Thoracic aortic aneurysm (TAA) | 5,704 | 460,490 | I9_THAORTANEUR |
| Aortic aneurysm overall (AA) | 10,694 | 460,490 | I9_AORTANEUR |
| Aortic dissection (AD) | 1,219 | 460,490 | I9_AORTDIS |

## Key results

- Seven taxon–outcome associations reached nominal significance in the discovery cohort (two for AAA, two for AA, three for AD) and six taxa reached nominal significance in AAAgen.
- No taxon survived Benjamini–Hochberg FDR correction in either cohort.
- No taxon replicated for AAA (0/140).
- Cross-cohort estimates correlated weakly (Pearson r = 0.20; 77/140 concordant direction).
- Nine of the 13 nominally significant cohort–outcome combinations lost significance when any single instrument was removed.
- Reverse MR yielded 34 estimates; none was significant (minimum P = 0.18).

## Repository contents

```text
.
├── analysis/
│   ├── mr_R13_AAAgen_dual_v2.R    # main forward/reverse MR and sensitivity analyses
│   └── make_figures.R             # Figures 1-3 and Supplementary Figures S1-S3
├── data/
│   ├── harmonised_dat.csv         # harmonised MiBioGen exposure instruments
│   ├── snps_hg19.csv              # hg19 coordinates for exposure SNPs
│   └── README.md                  # data sources and download instructions
├── results/
│   ├── forward_R13_ivw_fdr.csv
│   ├── forward_AAAgen_ivw_fdr.csv
│   ├── cross_cohort_lookup.csv
│   ├── sensitivity_all.csv
│   ├── loo_summary.csv
│   ├── presso_candidates.csv
│   ├── reverse_all.csv
│   ├── Additional_file_1_Tables_S1-S8.xlsx
│   └── README.md
├── CITATION.cff
├── LICENSE
├── .zenodo.json
└── ZENODO_UPLOAD_GUIDE.md
```

## Data sources

- **MiBioGen**: https://mibiogen.gcc.rug.nl
- **FinnGen release 13**: https://www.finngen.fi/en/access_results
- **FinnGen R13 manifest**: https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv
- **AAAgen**: https://csg.sph.umich.edu/willer/public/AAAgen2023/

Raw GWAS files are not included because of size. `data/harmonised_dat.csv` and `data/snps_hg19.csv` are included so the exposure side can be checked directly.

## How to reproduce

1. Install R 4.3.1 or later.
2. Install R packages:

```r
install.packages(c("data.table", "dplyr", "ggplot2", "cowplot"))
# TwoSampleMR and MRPRESSO are available from GitHub / CRAN; versions used:
# TwoSampleMR 0.4.26, MRPRESSO 1.0
```

3. Download the FinnGen R13 files for the four endpoints into `finngen-public-data-r13/`.
4. Download the AAAgen meta-analysis file into `AAAGen/`.
5. Put `harmonised_dat.csv` and `snps_hg19.csv` in the repository root, or copy them from `data/`.
6. Run:

```r
setwd("path/to/this/repository")
source("analysis/mr_R13_AAAgen_dual_v2.R")
```

7. Figures:

```r
setwd("path/to/this/repository")
source("analysis/make_figures.R")
```

## Important notes

- AAAgen is an **independent external replication cohort**, not a clinical validation cohort.
- The analysis uses summary-level data only; no individual-level data were accessed.
- The reverse MR analysis is exploratory because disease instruments were drawn from variants present in the exposure dataset.
- The study was approved by the Ethics Committee of Guilin Medical University (approval no. 2023-KY (0106)).

## License

Code is released under the MIT License (see `LICENSE`). The source GWAS data remain under the terms of the original data providers.
