# TCGA-LUAD Copper-Related Prognostic Signature

This repository contains selected analysis scripts, figures, and aggregate results associated with a copper-related prognostic signature study in TCGA-LUAD.

## Six-gene signature
RHOV, MT1H, LOXL2, SOD3, SLC31A2, and MT1A. The final stored analysis comprises 504 unique patients.

## Scripts
- `scripts/01_Reconstruct_6Gene_Risk_Score.R`: verifies the stored risk score from log2(x+1)-transformed expression and recovered coefficients. It **does not** reproduce the original LASSO training.
- `scripts/U5_Time_Dependent_ROC.R`: time-dependent ROC analysis. Check input paths before running.
- `reproducibility/U6_Software_Versions.R`: software version reporting.
- `scripts/archive/09_LASSO_Cox_Copper_Signature_HISTORICAL.R`: historical 13-candidate-gene Cox-LASSO workflow. **Not validated as the exact final 504-patient training procedure; do not cite as the final model training code.**

## Data availability and limitations
The score reconstruction script expects a local file at `drug_sensitivity/results/U2_6Gene_504/TCGA_504_Final_Risk_Groups.rds`. That input is **not included** in this public package. To run the script, obtain the appropriately processed data under applicable TCGA/GDC access and redistribution terms. The stored risk score can be reconstructed to numerical precision, but the original fitted LASSO model and full training provenance are not archived here. Some analyses are represented by output figures/tables rather than fully reproducible scripts.

## Results
Selected LASSO cross-validation figure, drug sensitivity summary, expression heatmap, time-dependent ROC figure and verified ROC summary. Results are research outputs, not independent external validation.

## Reproducibility
Run R scripts from the project root, install the packages referenced in the scripts, and inspect `reproducibility/` for recorded software versions. Patient-level files should not be publicly uploaded without a data-sharing review.
