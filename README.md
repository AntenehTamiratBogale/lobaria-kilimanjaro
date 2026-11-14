# lobaria-kilimanjaro

Analysis code and processed inputs for **A shared bacterial core and partner-specific functional variation in *Lobaria pulmonaria* across natural populations**.

This study combines 16S rRNA gene amplicon profiles and metaproteomics from 21 lichen thalli: three biological replicates in each of seven populations (Mweka1, Mweka2, Mweka3, Marangu, Maua, Umbwe and Styria). Styria is one designated reference population. The comparisons concern the sampled populations.

## What is included

- `Lobaria_complete_JPG_only.R`: complete R analysis and figure workflow.
- `input/`: four processed CSV files used by the analysis.
- `python_results/`: independent reference calculations and shared permutation/bootstrap draws from the earlier Python reanalysis. These files are reference inputs for cross-language checking, rather than a replacement for native R inference.
- `reference_results/native_R/`: selected previously completed native R results and verification tables.
- `scripts/`: dependency installation, input checks and explicit full-analysis/figure entry points.
- `docs/`: methods, figure index, data dictionary, validation status and provenance.

Raw sequencing reads and mass-spectrometry RAW files belong in appropriate data repositories and are not included in this code repository. No accession from another Lobaria study is assigned to this dataset here.

## Quick start

Use R 4.4 or later. The earlier native analysis was completed with R 4.4.1, vegan 2.6-8, limma 3.60.4 and DESeq2 1.44.0. The current package installer uses the packages available for your R/Bioconductor installation; it does not provide an exact frozen environment.

Open the repository root as your working directory, then run:

```r
source("scripts/install_packages.R")
source("scripts/check_inputs.R")
source("scripts/run_analysis.R")
```

From a terminal in the repository root:

```sh
Rscript scripts/install_packages.R
Rscript scripts/check_inputs.R
Rscript scripts/run_analysis.R
```

The full analysis can take substantial time. It fits native DESeq2 and limma models and runs permutation, bootstrap and sensitivity analyses.

After a completed local run, regenerate the figures with:

```r
source("scripts/make_figures.R")
```

Alternatively, `source("Lobaria_complete_JPG_only.R")` uses automatic mode: it reuses completed `R_results` objects when available and runs the full workflow when they are absent. The figure workflow also fits the added whole-metaproteome OG exclusion model needed for panel C of Supp-F8.

## Outputs

Results are created in `R_results/`. Manuscript graphics are exported as **600-dpi JPG only** in:

- `R_results/Final_Manuscript_JPG/Main/Figure1.jpg` through `Figure5.jpg`.
- `R_results/Final_Manuscript_JPG/Supplementary/Supp-F1.jpg` through `Supp-F21.jpg`.

Run-specific captions and numbering are written to `FIGURE_CAPTIONS.txt` and `FIGURE_INDEX.csv` in the same figure directory. Generated outputs are excluded from Git by default. The processed inputs and selected reference results are tracked.

See [the figure index](docs/FIGURES.md) and [analysis methods](docs/METHODS.md).

## Interpretation

The strict Kilimanjaro bacterial core consists of ASVs detected with at least one read in all 18 Kilimanjaro thalli. The previously completed analysis identified 118 core ASVs; 92 were detected in all three Styria thalli and 109 in at least one. The 109 describe broader recovery, not a strict shared core across all samples. Figures emphasize 118 and 92.

At OG resolution including Cyanobacteria, each Kilimanjaro population centroid was farther from Styria than from any other Kilimanjaro centroid in bacterial and whole-metaproteome analyses. The exact reference-relabelling test has seven choices and minimum attainable P = 1/7 = 0.143. Observed ordering and conventional outlier significance are separate conclusions.

Relative NSAF describes detected protein representation. Family-level protein annotations do not link proteins to exact ASVs. Nonsignificance does not establish equivalence. The selected Marangu-Styria volcano contrast is exploratory: Marangu is selected once by maximum all-bacterial OG centroid distance to Styria and retained across assays and exclusion sensitivity.

## Validation and manuscript status

The earlier returned native R run completed successfully. All 21 supplied Python-R checks passed; 77,234 limma coefficients matched independent calculations, with additional checks of BH adjustment and protein detection eligibility. See [validation status](docs/VALIDATION.md).

The latest revised graphics and the added whole-metaproteome exclusion model still require execution and visual review on the target R installation. This repository preparation does not certify new model results. The manuscript title, author declarations, collection-date discrepancies and dataset accessions must be finalized by the authors before a submission release.

## Citation and reuse

See `CITATION.cff` for repository citation metadata. The manuscript citation and complete author list can be added when finalized. A reuse license has not been selected; no license is inferred for third-party code or data. Contact the manuscript authors for permission until a license is explicitly added.
