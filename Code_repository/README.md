# Gut microbiota and aortic aneurysm/dissection: code and summary results

This repository accompanies the manuscript:

> **Gut microbiota and aortic aneurysm/dissection: a two-sample Mendelian randomization study with independent replication**

The analysis code, derived data, and summary result tables are in the [`Code_repository/`](Code_repository/) directory.

## What is in `Code_repository/`

```text
Code_repository/
├── analysis/
│   ├── mr_R13_AAAgen_dual_v2.R    # main forward/reverse MR and sensitivity analyses
│   └── make_figures.R             # Figures 1-3 and Supplementary Figures S1-S3
├── data/
│   ├── harmonised_dat.csv         # harmonised MiBioGen exposure instruments
│   ├── snps_hg19.csv              # hg19 coordinates for exposure SNPs
│   └── README.md
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
├── README.md
└── ZENODO_UPLOAD_GUIDE.md
```

## Study overview

| Item | Description |
|---|---|
| Exposure | MiBioGen gut microbiome GWAS (n = 18,340; 140 bacterial taxa; 1,363 unique SNPs) |
| Discovery | FinnGen release 13, four aortic outcomes |
| Replication | AAAgen abdominal aortic aneurysm GWAS (39,221 cases; 1,086,107 controls) |
| Primary method | Inverse-variance weighted MR |
| Software | R 4.3.1; TwoSampleMR 0.4.26; MRPRESSO 1.0 |

## Data sources

- MiBioGen: https://mibiogen.gcc.rug.nl
- FinnGen release 13: https://www.finngen.fi/en/access_results
- AAAgen: https://csg.sph.umich.edu/willer/public/AAAgen2023/

Raw GWAS files are not included because of size. The harmonised exposure data and hg19 coordinates are provided so the exposure side can be checked directly.

## Code DOI

`[Zenodo DOI to be inserted after deposit]`

## License

Code is released under the MIT License. The source GWAS data remain under the terms of the original data providers.
