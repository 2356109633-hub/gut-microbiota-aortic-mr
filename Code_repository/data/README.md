# Data inputs

## Included files

| File | Description |
|---|---|
| `harmonised_dat.csv` | Harmonised MiBioGen exposure instruments used by the main script |
| `snps_hg19.csv` | hg19/GRCh37 coordinates for exposure SNPs, used to match the AAAgen dataset |

## External files to download

### FinnGen release 13

Download the four endpoint files and place them in a folder named `finngen-public-data-r13/` in the repository root:

- `finngen_R13_I9_ABAORTANEUR.gz`
- `finngen_R13_I9_THAORTANEUR.gz`
- `finngen_R13_I9_AORTANEUR.gz`
- `finngen_R13_I9_AORTDIS.gz`

Source: https://www.finngen.fi/en/access_results  
Manifest: https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv

### AAAgen

Download the AAAgen meta-analysis summary statistics and place them in a folder named `AAAGen/` in the repository root. The script accepts a `.txt`, `.tsv` or `.gz` file inside that folder.

Source: https://csg.sph.umich.edu/willer/public/AAAgen2023/

## Notes

- AAAgen is reported on hg19 (GRCh37) without rsIDs. Variants are matched by position, assigned the corresponding exposure rsID and allele-aligned before harmonisation.
- Raw GWAS files are not redistributed here because of their size; please obtain them from the original providers.
- No individual-level data are included.
