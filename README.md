# Persistent postinfectious cough: reproducibility repository

Analysis code for *From biological signals to treatable traits in persistent postinfectious cough: a quantitative translational stress test*.

## Repository contents

| Path | Purpose |
|---|---|
| `01_proteomics_models.py` | Cohort preparation, circulating-factor models and interval compatibility |
| `02_prospective_power_tost.py` | Detection, exclusion, TOST equivalence and interval-width assurance |
| `03_airway_transcriptomics.R` | Five transcriptomic contexts, ROAST and exact permutation analyses |
| `04_laryngeal_and_trials.py` | Laryngeal physiology, participant responses and published trial contrasts |
| `05_generate_figures.py` | Figures 1–5 and Figure S1 |
| `data/raw/proteomics/` | Original INCOV and MY-LC workbooks |
| `data/airway/` | Transcript quantifications, count matrices, sample metadata and gene annotations |
| `data/clinical/` | Published laryngeal supplements, individual trial responses and aggregate trial data |
| `data/figure_metadata/` | Evidence-source relationships and measurement availability |
| `data/input_manifest.tsv` | Input paths, byte counts, MD5 and SHA-256 checksums |
| `data/input_sources.tsv` | Input definitions and source publications or accessions |
| `data/processed/proteomics/` | Derived cohort data and eligibility records |
| `output/tables/` | Numerical results |
| `output/figures/` | TIFF, PDF, SVG and preview images |
| `assets/fonts/` | Optional local installation of licensed Times New Roman fonts |

## Environment

The Python modules use Python 3.12 with NumPy 2.3.5, pandas 2.2.3, SciPy 1.17.0, openpyxl 3.1.5, Matplotlib 3.10.8 and Pillow 13.0.0. Module 03 uses R 4.4.3, edgeR 4.4.0, limma 3.62.1, tximport 1.34.0, statmod 1.5.2 and jsonlite 2.0.0. ROAST is executed through its native Bioconductor implementation.

```bash
python -m venv .venv
.venv/bin/python -m pip install numpy==2.3.5 pandas==2.2.3 scipy==1.17.0 openpyxl==3.1.5 matplotlib==3.10.8 Pillow==13.0.0
```

On Windows, the Python executable is `.venv/Scripts/python.exe`.

A compatible R environment can be installed with micromamba:

```bash
micromamba create -y -p .r-env -c conda-forge -c bioconda r-base=4.4.3 bioconductor-edger=4.4.0 bioconductor-limma=3.62.1 bioconductor-tximport=1.34.0 r-statmod=1.5.2 r-jsonlite=2.0.0
```

Times New Roman must be installed on the rendering system, or its licensed `.ttf` files placed in `assets/fonts/`. Fonts are not distributed. Figure export stops if the requested font is unavailable.

## Data acquisition

The accompanying `Reproduction_Input_Data.tar.gz` contains the input files under their repository-relative paths. Extract it at the repository root:

```bash
python -m tarfile -e Reproduction_Input_Data.tar.gz .
```

The input manifest records the exact files used. Module 01 can also retrieve the three original proteomic workbooks with `--download`. Each source retains its original attribution and reuse conditions.

| Analysis | Public source | Analysis input |
|---|---|---|
| INCOV | [Talla et al.](https://doi.org/10.1038/s41467-023-38682-4); [source repository](https://github.com/aifimmunology/PASC-proteomics-talla-vasaikar-et-al/tree/e34c6de8583983074251e320acc7fd904c009739/Input%20data) | `JimHeath_proteomics_metadata.xlsx`; `JimHeath_proteomics.xlsx` |
| MY-LC | [Klein et al.](https://doi.org/10.1038/s41586-023-06651-y) | `41586_2023_6651_MOESM4_ESM.xlsx` |
| Paired epithelial culture | [GSE286429](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286429) | Twelve Salmon quantifications, GENCODE v44 transcript-to-gene mapping and donor manifest |
| Longitudinal nasal epithelium | [GSE261777](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261777) | Published count matrix and metadata identifying six paired donors |
| Radiological comparison | [GSE261777](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261777) | The same count matrix, published clinical labels and deconvolution fractions |
| Nasal donor pseudobulks | [GSE232623](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232623) / [GSE232627](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232627) | Human gene counts summed within biological donors, cell-quality measurements and donor covariates |
| Whole-blood cough annotation | [GSE226260](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226260) | Counts and clinical annotation for one indexed post-acute sample per participant |
| Laryngeal physiology | [García-Vicente et al.](https://doi.org/10.1371/journal.pone.0283758) | Original supplementary files `pone.0283758.s002.docx`, `s003.docx` and `s004.xlsx` |
| Amitriptyline trial | [Jeyakumar et al.](https://doi.org/10.1097/01.mlg.0000244377.60334.e3) | All 28 individual responses published in Table III |
| Endpoint benchmarks | Source citations in `data/input_sources.tsv` and the clinical extraction tables | Published between-arm effects, within-arm changes and endpoint classifications |

The transcriptomic entry points are Salmon quantifications and gene-count matrices. The distributed inputs include the donor-level pseudobulks used for modelling; the five analysis modules do not repeat FASTQ processing or cell-level read alignment. Source metadata document the figure-derived radiological labels, donor-covariate linkage, and whole-blood symptom annotations.

## Execution

Run all commands from the repository root:

```bash
.venv/bin/python 01_proteomics_models.py --download
.venv/bin/python 02_prospective_power_tost.py --validate
micromamba run -p .r-env Rscript 03_airway_transcriptomics.R
.venv/bin/python 04_laryngeal_and_trials.py
.venv/bin/python 05_generate_figures.py
```

Omit `--download` when the proteomic workbooks are present. Module 02 is independent of participant-level data. Module 04 is independent of modules 01–03. Module 05 reads their numerical outputs.

Module 03 accepts `culture`, `longitudinal`, `radiology`, `pseudobulk`, `blood` or `summaries` as optional positional arguments. Running without arguments executes all five contexts and the summary/permutation stage. The radiology module uses the longitudinal gene-statistics output. `summaries` requires the completed context outputs. Input mismatches or numerical errors terminate execution.

## Statistical specification

**Circulating factors.** Cohort eligibility uses documented cough, assessment time and the specified covariates. Released assay values are standardised within each model population; zero and negative scores are retained. Ordinary least squares uses HC3 covariance and a Student t reference with residual degrees of freedom. Benjamini–Hochberg correction covers twelve primary comparisons, twelve comparisons within each sensitivity specification, and separate exploratory families of 448 INCOV and 138 MY-LC assays. The 586 exploratory assays include nonprotein soluble factors.

The four sensitivity specifications use rank-normalised outcomes, assessment at least 56 days after symptom onset, cohort-specific respiratory-history restrictions and additional dyspnoea adjustment. Each participant omission repeats outcome standardisation and model fitting. Compatibility is evaluated against positive and negative differences of 0.2, 0.5 and 0.8 assay standard deviations using nominal and simultaneous confidence intervals. The expected outputs comprise twelve primary, 48 sensitivity and 1,212 participant-omission models.

**Prospective designs.** Independent normal groups are assumed to have equal unknown variance. Detection uses the noncentral t distribution. Exclusion evaluates a positive margin under a zero true difference, using either a one-sided upper bound or the upper endpoint of a two-sided interval. TOST power uses normal–chi-square integration; interval-width assurance is calculated separately. Integer group sizes use `ceiling(p*N)` and `N-ceiling(p*N)`. Recruitment allowances apply 20% attrition separately to each group. Minimum totals are checked against the preceding integer total. These calculations concern independent-group designs; covariate-adjusted HC3 models have a different sampling structure.

**Transcriptomics.** Gene filtering uses `filterByExpr`; TMM normalisation precedes voom, robust limma moderation and 19,999-rotation ROAST tests. Directional and mixed tests form separate ten-programme families. Gene identifiers, representative-gene selection and coverage rules are retained for each context. Culture counts are imported with `countsFromAbundance="lengthScaledTPM"`. Culture models estimate the between-group difference in within-donor stimulation response; longitudinal models estimate late-minus-early changes within donor. Radiological sensitivity models include two cell-composition principal components. Pseudobulk models use biological donors, with compartment-specific and cell-quality sensitivities. Blood models include age, sex and elapsed time.

Donor omissions retain the specified gene universe and repeat normalisation and model fitting. The longitudinal score analysis enumerates all 64 donor-level sign patterns. Unweighted studentised score analyses enumerate 2,324,784 radiological, 462 pseudobulk and 6,188 blood-label arrangements. The pseudobulk cell-composition comparison uses a separate 22-test family. RNG settings, seeds, design matrices, gene mappings and software versions are recorded in module 03 and its outputs.

**Laryngeal physiology and trial endpoints.** The 38 clinical records are linked to the published LEMG descriptions by participant identifier. Wilson intervals retain all participants in each specified subgroup. Five exploratory symptom comparisons use two-sided Fisher tests and one Benjamini–Hochberg family. Trial responses of at least 50% and 75% are reported separately. Risk differences use uncorrected Newcombe–Wilson intervals, with the zero control-response cell retained. Published adjusted LCQ and gefapixant effects are extracted without pooling. Procaterol within-arm intervals are calculated from reported means and standard deviations using a t reference; they are separate from its published repeated-measures treatment coefficient.

**Figures.** Figures 1–5 contain six panels each; Figure S1 contains two. Export formats are 600-dpi LZW-compressed TIFF, PDF, editable-text SVG and PNG previews. Numerical labels are rounded only during presentation.

## Outputs

| Location | Contents | Manuscript correspondence |
|---|---|---|
| `data/processed/proteomics/` | Clinical data, eligibility and assay mapping | Cohort definitions and Table S11 |
| `output/tables/proteomics/primary_panel_results.csv` | Twelve primary comparisons | Table S2; Figure 2 |
| `output/tables/proteomics/sensitivity_results.csv` | Forty-eight sensitivity models | Table S3 |
| `output/tables/proteomics/leave_one_out_results.csv` and `loo_summary.csv` | Participant-omission estimates and summaries | Tables S4a–S4b |
| `output/tables/proteomics/exploratory_assay_results.csv` | 586 exploratory comparisons | Table S5; Figure 2 |
| `output/tables/proteomics/primary_precision.csv` and `effect_bound_compatibility.csv` | Interval precision and compatibility | Tables S6a–S6b; Figure 2 |
| `output/tables/proteomics/Table1_numeric_long.csv` | Model-population characteristics | Table 1 |
| `output/tables/prospective/prospective_design_all_scenarios.csv` | 192 detection, exclusion and equivalence scenarios | Design calculations |
| `output/tables/prospective/Table_S13_selected_scenarios.csv` | Six selected design scenarios | Table S13 and Table 2 |
| `output/tables/prospective/prospective_CI_width_assurance.csv` | 36 interval-width scenarios | Table S14 |
| `output/tables/prospective/N100_detection_power.csv` | Four N=100 scenarios | Table S15 |
| `output/tables/prospective/prospective_design_curves.csv` | 1,888 design probabilities | Figure S1 |
| `output/tables/airway/primary_programmes.tsv` | Fifty primary programme comparisons | Tables S7 and S12; Figure 3 |
| `output/tables/airway/culture/` | Culture interaction, secondary contrasts, sensitivities and donor omissions | Table S12 |
| `output/tables/airway/longitudinal/` | Paired ROAST, donor changes, 64-pattern sign-flip results and omissions | Tables S7 and S12; Figure 3 |
| `output/tables/airway/radiology/` | Radiological ROAST, composition adjustments, control omissions and exact score tests | Tables S7 and S12; Figure 3 |
| `output/tables/airway/pseudobulk/` | Donor and compartment models, omissions, score and cell-composition tests | Table S12; Figure 3 |
| `output/tables/airway/blood/` | Cough-annotation models, sensitivities, omissions and score tests | Table S12 |
| `output/tables/clinical/clinical_lemg_merged.csv` | Linked clinical and electrophysiological records | Tables S9a–S9c; Figure 4 |
| `output/tables/clinical/lemg_descriptive_proportions.csv` | Subgroup proportions and Wilson intervals | Table S9a; Figure 4 |
| `output/tables/clinical/exploratory_within_case_symptom_contrasts.csv` | Five Fisher tests and adjusted probabilities | Table S9b; Figure 4 |
| `output/tables/clinical/trial_patient_response.csv` and `trial_reanalysis.json` | Individual responses and both response thresholds | Table S10a; Figure 5 |
| `output/tables/clinical/between_arm_effects.csv` and `within_arm_changes.csv` | Adjusted effects and arm-level changes | Tables S10b–S10c; Figure 5 |
| `output/tables/clinical/gefapixant_effects.csv` and `trial_endpoint_matrix.csv` | Published relative effects and endpoint classifications | Figure 5 |
| `output/figures/Figure1`–`Figure5`, `FigureS1` | Six figure sets | Main and supplementary figures |

Each analysis records its software environment. Module 02 with `--validate` writes Monte Carlo comparisons in `simulation_checks.csv`. Numerical tables retain full computational precision. The figure manifest records file hashes and the font identity.

## Methodological references

- MacKinnon JG, White H. Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties. *J Econometrics*. 1985;29:305–325. doi:10.1016/0304-4076(85)90158-7.
- Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B*. 1995;57:289–300. doi:10.1111/j.2517-6161.1995.tb02031.x.
- Shieh G. Exact power and sample size calculations for the two one-sided tests of equivalence. *PLoS One*. 2016;11:e0162093. doi:10.1371/journal.pone.0162093.
- Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*. 2010;26:139–140. doi:10.1093/bioinformatics/btp616.
- Law CW, Chen Y, Shi W, et al. voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biol*. 2014;15:R29. doi:10.1186/gb-2014-15-2-r29.
- Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res*. 2015;43:e47. doi:10.1093/nar/gkv007.
- Wu D, Lim E, Vaillant F, et al. ROAST: rotation gene set tests for complex microarray experiments. *Bioinformatics*. 2010;26:2176–2182. doi:10.1093/bioinformatics/btq401.
- Wilson EB. Probable inference, the law of succession, and statistical inference. *J Am Stat Assoc*. 1927;22:209–212. doi:10.1080/01621459.1927.10502953.
- Newcombe RG. Interval estimation for the difference between independent proportions: comparison of eleven methods. *Stat Med*. 1998;17:873–890. doi:10.1002/(sici)1097-0258(19980430)17:8<873::aid-sim779>3.0.co;2-i.
