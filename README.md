# MSSA precipitation regionalization of North Sumatra

Code and data for "Precipitation regionalization of North Sumatra using multichannel singular spectrum analysis" (ICDSOS 2026).

| File | Purpose |
|---|---|
| `gee_chirps_sumut_v2.js` | Google Earth Engine: daily CHIRPS areal means per regency/city, 1995–2024 |
| `Precip_CHIRPS_Sumut_33kab_1995_2024.csv` | Extracted data (10,958 days × 33 units, mm/day) |
| `zonasi_mssa_sumut.R` | Main analysis in R (Rssa): MSSA, W-correlation, Ward clustering, baselines |
| `zonasi_mssa_sumut.py` | Independent Python implementation of the main analysis |
| `rssa_like.py` | Python MSSA (horizontal trajectory matrix, as in Rssa) |
| `a1.py`–`a5.py` | Revision analyses: component selection per configuration, harmonic baseline, 122-day pair, transforms, robustness (15 configurations), bootstrap Jaccard, Oldeman and BMKG criteria, pixel counts |
| `figs_rev.py` | Figures 3–6 of the revised paper |
| `gee_peta_zona_sumut.js` | Google Earth Engine map of the regions |

Paths inside the scripts point to the authors' working folders; edit them before running.
