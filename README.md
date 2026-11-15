# Post-COVID symptom-specific biomarker analyses

Code and numerical results for *Precision, Phenotype Misclassification, and Preprocessing in Post–COVID-19 Cough Biomarker Studies*.

Version: 2.0.0.

## Contents

| Module | Analysis |
|---|---|
| `01_proteomics_models.py` | INCOV and MY-LC eligibility; six-protein OLS–HC3 models; sensitivity and participant-omission models; exploratory assays; interval compatibility |
| `02_prospective_power_tost.py` | Independent-group detection, exclusion, exact TOST and interval-width assurance |
| `03_airway_transcriptomics.R` | Five donor-level transcriptomic contrasts, ROAST, donor omissions and exact score tests |
| `04_ademp_simulation.py` | Empirical calibration, ADEMP simulation, implementation validation and conditional sample-size estimation |
| `05_generate_figures.py` | Main Figures 1–3 and supplementary eFigures 1–2 |

`data/raw/proteomics/` contains the released source workbooks. `data/airway/` contains the transcript quantifications, count matrices, donor metadata and gene annotations used for modelling. `data/ademp_calibration.json` contains the empirical simulation parameters. `data/Table_S1_phenotype_availability.csv` records the eleven-source phenotype assessment. Source attribution and input checksums are in `data/input_sources.tsv` and `data/input_manifest.tsv`.

The distributed entry points for transcriptomics are Salmon quantifications and donor-level matrices. FASTQ alignment and cell-level read processing are upstream of these modules. The five contrasts include four airway contexts and one whole-blood context; the two GSE261777 contrasts share a source cohort.

## Environment

The analyses use Python 3.12 with NumPy, pandas, SciPy, openpyxl, Matplotlib, Pillow and threadpoolctl. Tested versions are recorded in `requirements.txt` and the numerical execution records. R analyses use R 4.4.3, edgeR 4.4.0, limma 3.62.1, tximport 1.34.0, statmod 1.5.2 and jsonlite 2.0.0. The R script checks package availability and records `sessionInfo()`. Archived environment reports retain library names and versions without machine-specific directory prefixes.

```bash
python -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
micromamba create -y -p .r-env -c conda-forge -c bioconda r-base=4.4.3 bioconductor-edger=4.4.0 bioconductor-limma=3.62.1 bioconductor-tximport=1.34.0 r-statmod=1.5.2 r-jsonlite=2.0.0
```

On Windows, the Python executable is `.venv/Scripts/python.exe`. The multiprocessing lock uses POSIX `fcntl`; simulation execution requires Linux, macOS or Windows Subsystem for Linux. Times New Roman must be installed for figure export; licensed font files can also be placed in `assets/fonts/`. Fonts are not redistributed.

## Public data

| Input | Source |
|---|---|
| INCOV clinical and proteomic workbooks | Talla et al., doi:10.1038/s41467-023-38682-4; https://github.com/aifimmunology/PASC-proteomics-talla-vasaikar-et-al |
| MY-LC clinical and circulating factors | Klein et al., Supplementary Table 3; doi:10.1038/s41586-023-06651-y |
| Paired epithelial cultures | GEO GSE286429 |
| Longitudinal nasal and radiological contrasts | GEO GSE261777 |
| Nasal epithelial donor pseudobulks | GEO GSE232623 and GSE232627 |
| Whole-blood symptom annotation | GEO GSE226260 |

Input files retain their source attribution and reuse conditions. The phenotype inventory describes the inspected public releases, not all information held by the original investigators. Missing symptom listings are not treated as confirmed negative responses.

## Execution

Run commands from the repository root. The first three modules produce the empirical outputs:

```bash
.venv/bin/python 01_proteomics_models.py
.venv/bin/python 02_prospective_power_tost.py --validate
micromamba run -p .r-env Rscript 03_airway_transcriptomics.R
```

Module 01 accepts `--download` to retrieve absent proteomic workbooks. Module 03 accepts `culture`, `longitudinal`, `radiology`, `pseudobulk`, `blood` or `summaries`; running without an argument performs the complete sequence. The radiological analysis depends on the longitudinal gene-statistics output, and `summaries` requires the completed context outputs.

The supplied calibration permits direct reproduction of the simulations. `prepare` can regenerate calibration after module 01:

```bash
.venv/bin/python 04_ademp_simulation.py prepare --source .
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python 04_ademp_simulation.py validate
.venv/bin/python 04_ademp_simulation.py run --workers 6
.venv/bin/python 04_ademp_simulation.py design --workers 6
.venv/bin/python 04_ademp_simulation.py summarise
.venv/bin/python 05_generate_figures.py
```

`validate` must precede `run` and must match the code and calibration. Completed scenario checkpoints are reused only with the same source, calibration and repetition settings. Parallel scheduling does not alter the random-number streams. Worker count can be reduced without changing the simulations. The full run uses 312 combinations with 10,000 repetitions and 96 additional design combinations with 20,000 repetitions, totalling 5,040,000 target datasets. Auxiliary source samples and validation datasets are separate.

All aggregate results are distributed. Replicate-level checkpoint arrays are generated by the simulation commands and omitted from the publication archive. Running `summarise` requires those checkpoints. Figures can be generated directly from the supplied aggregate tables without repeating the simulations.

## Statistical specification

### Empirical analyses

The primary proteomic population comprises 105 INCOV participants (26 recorded cough) and 97 MY-LC participants (24 recorded cough). Released assay values are standardised within complete-model observations. Models adjust for cough, age, sex, BMI, log elapsed time and acute hospitalisation. HC3 covariance uses a Student-t reference with residual degrees of freedom. BH correction covers twelve primary comparisons, twelve comparisons within each sensitivity specification, and separate 448-assay and 138-assay exploratory panels. Four sensitivity specifications yield 48 models; participant omission yields 1,212 models. Nominal and Bonferroni intervals are compared with reference values of 0.2, 0.5 and 0.8 assay SD.

Transcriptomic analyses use `filterByExpr`, TMM, voom–limma and 19,999-rotation ROAST. Directional and mixed tests have separate ten-programme correction families. Longitudinal exact sign flips enumerate 64 patterns. Studentised score tests enumerate 2,324,784 radiological, 462 pseudobulk and 6,188 blood-label arrangements. Unweighted score tests and precision-weighted ROAST summaries are retained separately.

### Analytical designs

Independent normal groups have equal unknown variance. Detection uses the noncentral t distribution; exact TOST integrates over independent normal and chi-square variables. Exclusion and interval-width assurance are separate objectives. Integer group sizes use `ceiling(p*N)`; recruitment allowances apply 20% attrition within each group. The 466, 169 and 509 totals are analytical references, not general requirements for the adjusted models.

### ADEMP simulation

The main generator resamples joint covariate vectors within the two empirical cough groups and retains group-specific six-protein residual covariance. Latent prevalence is 25%; sensitivity/specificity affects only the observed phenotype. Prespecified pairs are 1.00/1.00, 0.80/1.00, 1.00/0.95, 0.80/0.95 and 0.80/0.85. Sample sizes are 100, 169, 466, 509 and 1000. Effect vectors are all null, IFNG alone at +0.30 SD, six +0.30-SD effects or six +0.50-SD effects. Seven extensions assess variance equality, nuisance R², residual distributions, source interactions, prevalence and negative boundaries.

Methods comprise an oracle model, joint adjustment, same-sample residualisation with readjustment, outcome-only residualisation followed by group comparison, and source residualisation. Native and subgroup-standardised coefficients are both retained. Effect targets and margins are translated with the realised SD when scales change. The reference SD is the pooled within-group outcome SD at zero phenotype effect and latent prevalence 25%.

Outputs distinguish coverage of the imposed latent effect from coverage of the observed-phenotype projection parameter. Exclusion requires the upper 95% interval limit below a positive margin; TOST requires the 90% interval within symmetric margins. BH FDR is evaluated with all-null and sparse alternatives. Failed analyses remain in decision-rate denominators with no successful decision; coefficient and coverage summaries condition on valid fits.

Every probability and estimator summary includes a Monte Carlo standard error. Six-protein mean rates are averaged within each repetition before estimating their Monte Carlo uncertainty. Conditional N80 uses monotone interpolation rounded upwards to ten participants, with pointwise Monte Carlo brackets retained. These brackets are not simultaneous confidence intervals for the adaptively selected sample size. Misclassified-phenotype N80 describes detection of a proxy association, not restoration of latent-effect coverage.

## Output correspondence

| Location | Manuscript item |
|---|---|
| `data/Table_S1_phenotype_availability.csv` | eTable 1; source selection |
| `output/tables/proteomics/Table1_numeric_long.csv` | Table 1 |
| `output/tables/proteomics/primary_panel_results.csv` | Figure 1; eTable 2 |
| `output/tables/proteomics/sensitivity_results.csv` | eTable 3 |
| `output/tables/proteomics/leave_one_out_results.csv`, `loo_summary.csv` | eTable 4 |
| `output/tables/proteomics/exploratory_assay_results.csv` | Figure 1; eTable 5 |
| `output/tables/proteomics/primary_precision.csv`, `effect_bound_compatibility.csv` | Figure 1; eTable 6 |
| `output/tables/airway/primary_programmes.tsv` | Figure 2; eTable 7 |
| `output/tables/airway/longitudinal/` | Figure 2; eTable 8 |
| `output/tables/airway/` context directories | eTable 9 and complete transcriptomic sensitivities |
| `data/processed/proteomics/primary_assay_map.csv` | eTable 10 |
| `output/tables/prospective/selected_design_scenarios.csv`, `prospective_design_all_scenarios.csv` | eTable 11; analytical part of Table 2 |
| `output/tables/prospective/prospective_CI_width_assurance.csv`, `N100_detection_power.csv` | eTable 12 |
| `output/ademp/empirical_calibration.csv`, `*_residual_correlation.csv` | Figure 3A; eTable 13 |
| `output/ademp/design.json`, `output/design/conditional_design_summary.json` | eTable 14 |
| `output/ademp/simulation_estimates.csv`, `simulation_panels.csv`, `simulation_bh12.csv` | Figure 3; eTables 15–17 |
| `output/ademp/key_findings.csv` | Main simulation results; eTable 16 |
| `output/ademp/paired_method_comparisons.csv`, `misclassification_projection.csv` | eTable 17 |
| `output/design/Table_2_conditional_design.csv`, `Table_2_design_ranges.csv` | Table 2; Figure 3F; eTable 18 |
| `output/ademp/validation.json`, `computational_checks.json` | eTable 19 |
| `output/figures/Figure_1`–`Figure_3` | Main figures |
| `output/figures/Figure_S1`, `Figure_S2` | Supplementary figures |
| `output/supplementary_table_index.csv` | Complete supplementary table-to-file index |

Figure export includes PDF, SVG, 600-dpi LZW TIFF and PNG previews. `SHA256SUMS.txt` covers the files distributed in this release, excluding the checksum file itself. Simulation checkpoint arrays and regenerable expression intermediates are not required to inspect the aggregate results.

## Methods references

- Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate statistical methods. *Stat Med*. 2019;38:2074–2102. doi:10.1002/sim.8086.
- MacKinnon JG, White H. Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties. *J Econometrics*. 1985;29:305–325. doi:10.1016/0304-4076(85)90158-7.
- Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B*. 1995;57:289–300. doi:10.1111/j.2517-6161.1995.tb02031.x.
- Shieh G. Exact power and sample size calculations for the two one-sided tests of equivalence. *PLoS One*. 2016;11:e0162093. doi:10.1371/journal.pone.0162093.
- Law CW, Chen Y, Shi W, et al. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biol*. 2014;15:R29. doi:10.1186/gb-2014-15-2-r29.
- Wu D, Lim E, Vaillant F, et al. ROAST: rotation gene set tests for complex microarray experiments. *Bioinformatics*. 2010;26:2176–2182. doi:10.1093/bioinformatics/btq401.
