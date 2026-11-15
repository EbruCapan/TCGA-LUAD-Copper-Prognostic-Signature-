# TCGA-LUAD: reconstruction of the stored six-gene risk score
# This script verifies a previously calculated score; it does not retrain LASSO-Cox.
# Run from the LUAD_Copper_Project root directory.

input_file <- file.path(
  "drug_sensitivity", "results", "U2_6Gene_504",
  "TCGA_504_Final_Risk_Groups.rds"
)
if (!file.exists(input_file)) {
  stop("Input not found: ", input_file,
       "\nSet the working directory to the LUAD_Copper_Project root.")
}

dat <- readRDS(input_file)
genes <- c("RHOV", "MT1H", "LOXL2", "SOD3", "SLC31A2", "MT1A")
required <- c("submitter_id", "risk_score", genes)
stopifnot(all(required %in% names(dat)))
stopifnot(nrow(dat) == 504L, !anyDuplicated(dat$submitter_id))

X <- as.matrix(dat[, genes, drop = FALSE])
stopifnot(is.numeric(X), all(is.finite(X)), all(X >= 0))
stopifnot(is.numeric(dat$risk_score), all(is.finite(dat$risk_score)))
X_log2 <- log2(X + 1)

# Coefficients reconstructed from the archived scores (not original glmnet object).
beta <- c(
  RHOV = 0.0950696240643004,
  MT1H = -0.0270633080447723,
  LOXL2 = 0.1395063363276147,
  SOD3 = -0.1255297650691920,
  SLC31A2 = -0.0950283370219768,
  MT1A = 0.0524027831990282
)
intercept <- -0.8105247919020858

reconstructed <- as.vector(intercept + X_log2[, names(beta), drop = FALSE] %*% beta)
difference <- reconstructed - dat$risk_score
max_diff <- max(abs(difference))
mean_diff <- mean(abs(difference))

cat("Patients:", nrow(dat), "\n")
cat("Maximum absolute difference:", format(max_diff, scientific = TRUE), "\n")
cat("Mean absolute difference:", format(mean_diff, scientific = TRUE), "\n")
stopifnot(max_diff < 1e-10)

output_dir <- file.path("results", "risk_score_validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
validation <- data.frame(
  Metric = c("Patients", "Maximum_absolute_difference", "Mean_absolute_difference"),
  Value = c(nrow(dat), max_diff, mean_diff)
)
write.csv(validation, file.path(output_dir, "validation_summary.csv"), row.names = FALSE)
cat("Risk score reconstruction verified.\n")
