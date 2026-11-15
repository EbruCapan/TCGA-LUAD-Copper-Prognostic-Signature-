library(survival)
library(glmnet)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(survminer)
library(timeROC)

counts <- readRDS("data_processed/TCGA_LUAD_primary_normal_counts.rds")
sample_info <- readRDS("data_processed/TCGA_LUAD_sample_info_filtered.rds")
clinical <- readRDS("clinical/TCGA_LUAD_clinical.rds")
deg_all <- openxlsx::read.xlsx("results/02_DEG_Tumor_vs_Normal_all.xlsx")

candidate_genes <- c(
  "RHOV", "CDKN2A", "CP", "STEAP1", "MT1H",
  "LOXL2", "STEAP2", "COA6", "MT1X",
  "SOD3", "SLC31A2", "MT1A", "MT1M"
)

sample_info$sample_id <- rownames(sample_info)

tumor_info <- sample_info %>%
  dplyr::filter(sample_type == "Primary Tumor")

tumor_counts <- counts[, tumor_info$sample_id]

# Ensembl ID temizleme
deg_all$gene_id_clean <- sub("\\..*", "", deg_all$gene_id)
count_gene_id_clean <- sub("\\..*", "", rownames(tumor_counts))

gene_map <- deg_all %>%
  dplyr::filter(gene_name %in% candidate_genes) %>%
  dplyr::select(gene_id_clean, gene_name) %>%
  dplyr::distinct()

idx <- match(gene_map$gene_id_clean, count_gene_id_clean)

gene_map <- gene_map[!is.na(idx), ]
idx <- idx[!is.na(idx)]

candidate_counts <- tumor_counts[idx, , drop = FALSE]
rownames(candidate_counts) <- gene_map$gene_name

# Aynı gen birden fazla satırdaysa ortalama al
candidate_counts_avg <- rowsum(candidate_counts, group = rownames(candidate_counts)) /
  as.vector(table(rownames(candidate_counts)))

expr_mat <- as.data.frame(t(candidate_counts_avg))
expr_mat$sample_id <- rownames(expr_mat)

existing_genes <- candidate_genes[candidate_genes %in% colnames(expr_mat)]

print(existing_genes)
print(length(existing_genes))

if (length(existing_genes) < 2) {
  stop("LASSO icin yeterli gen bulunamadi.")
}

clinical_surv <- clinical %>%
  dplyr::select(
    submitter_id,
    vital_status,
    days_to_death,
    days_to_last_follow_up
  )

clinical_surv$OS.time <- ifelse(
  is.na(clinical_surv$days_to_death),
  clinical_surv$days_to_last_follow_up,
  clinical_surv$days_to_death
)

clinical_surv$OS <- ifelse(
  clinical_surv$vital_status == "Dead",
  1,
  0
)

tumor_info$submitter_id <- substr(tumor_info$barcode, 1, 12)

model_data <- tumor_info %>%
  dplyr::select(sample_id, submitter_id) %>%
  dplyr::left_join(expr_mat, by = "sample_id") %>%
  dplyr::left_join(clinical_surv, by = "submitter_id") %>%
  dplyr::filter(!is.na(OS.time), OS.time > 0)

# 6. LASSO input
x <- as.matrix(model_data[, existing_genes])
x <- scale(log2(x + 1))

y <- Surv(model_data$OS.time, model_data$OS)

set.seed(123)

cvfit <- cv.glmnet(
  x,
  y,
  family = "cox",
  alpha = 1,
  nfolds = 10
)

best_lambda <- cvfit$lambda.min

lasso_model <- glmnet(
  x,
  y,
  family = "cox",
  alpha = 1,
  lambda = best_lambda
)

coef_mat <- as.matrix(coef(lasso_model))

lasso_genes <- data.frame(
  Gene = rownames(coef_mat),
  Coefficient = as.numeric(coef_mat[, 1])
) %>%
  dplyr::filter(Coefficient != 0)

# 7. Risk score
risk_score <- as.vector(
  x[, lasso_genes$Gene, drop = FALSE] %*% lasso_genes$Coefficient
)

model_data$risk_score <- risk_score

cutoff <- median(model_data$risk_score, na.rm = TRUE)

model_data$risk_group <- ifelse(
  model_data$risk_score >= cutoff,
  "High-risk",
  "Low-risk"
)

model_data$risk_group <- factor(
  model_data$risk_group,
  levels = c("Low-risk", "High-risk")
)

# 8. Excel ciktilari
openxlsx::write.xlsx(
  lasso_genes,
  "results/06_LASSO_Selected_Copper_Genes.xlsx",
  rowNames = FALSE
)

openxlsx::write.xlsx(
  model_data,
  "results/06_Copper_Risk_Score_TCGA.xlsx",
  rowNames = FALSE
)

# 9. LASSO CV curve
pdf("figures/06_LASSO_CV_Curve.pdf", width = 7, height = 6)
plot(cvfit)
dev.off()

# 10. Kaplan-Meier
fit <- survfit(
  Surv(OS.time, OS) ~ risk_group,
  data = model_data
)

p_km <- ggsurvplot(
  fit,
  data = model_data,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("#4DBBD5", "#E64B35"),
  title = "Copper-Related Risk Signature in TCGA-LUAD",
  legend.title = "Risk group",
  legend.labs = c("Low-risk", "High-risk")
)

ggsave(
  "figures/06_LASSO_RiskScore_KM.pdf",
  p_km$plot,
  width = 7,
  height = 5
)

ggsave(
  "figures/06_LASSO_RiskScore_KM.png",
  p_km$plot,
  width = 7,
  height = 5,
  dpi = 600
)

# 11. Time-dependent ROC
roc_res <- timeROC(
  T = model_data$OS.time,
  delta = model_data$OS,
  marker = model_data$risk_score,
  cause = 1,
  times = c(365, 1095, 1825),
  iid = TRUE
)

pdf("figures/06_TimeDependent_ROC.pdf", width = 7, height = 6)
plot(roc_res, time = 365, col = "red", title = FALSE)
plot(roc_res, time = 1095, add = TRUE, col = "blue")
plot(roc_res, time = 1825, add = TRUE, col = "darkgreen")
legend(
  "bottomright",
  legend = c(
    paste0("1-year AUC = ", round(roc_res$AUC[1], 3)),
    paste0("3-year AUC = ", round(roc_res$AUC[2], 3)),
    paste0("5-year AUC = ", round(roc_res$AUC[3], 3))
  ),
  col = c("red", "blue", "darkgreen"),
  lwd = 2
)
dev.off()

# 12. Risk score plot
risk_plot_df <- model_data %>%
  dplyr::arrange(risk_score) %>%
  dplyr::mutate(index = dplyr::row_number())

p_risk <- ggplot(
  risk_plot_df,
  aes(x = index, y = risk_score, color = risk_group)
) +
  geom_point(size = 1.8) +
  theme_bw() +
  labs(
    title = "Risk Score Distribution",
    x = "Patients",
    y = "Risk score"
  )

ggsave(
  "figures/06_RiskScore_Distribution.pdf",
  p_risk,
  width = 7,
  height = 5
)

# 13. Survival status plot
p_status <- ggplot(
  risk_plot_df,
  aes(x = index, y = OS.time, color = factor(OS))
) +
  geom_point(size = 1.8) +
  theme_bw() +
  labs(
    title = "Survival Status Distribution",
    x = "Patients",
    y = "Overall survival time",
    color = "Status"
  )

ggsave(
  "figures/06_Survival_Status_Distribution.pdf",
  p_status,
  width = 7,
  height = 5
)

# 14. Model performansi
cox_risk <- coxph(
  Surv(OS.time, OS) ~ risk_score,
  data = model_data
)

cox_summary <- summary(cox_risk)

summary_table <- data.frame(
  Model = "Copper Risk Score",
  HR = cox_summary$coefficients[, "exp(coef)"],
  lower95 = cox_summary$conf.int[, "lower .95"],
  upper95 = cox_summary$conf.int[, "upper .95"],
  p.value = cox_summary$coefficients[, "Pr(>|z|)"],
  AUC_1year = roc_res$AUC[1],
  AUC_3year = roc_res$AUC[2],
  AUC_5year = roc_res$AUC[3]
)

openxlsx::write.xlsx(
  summary_table,
  "results/06_LASSO_Model_Performance.xlsx",
  rowNames = FALSE
)

print(lasso_genes)
print(summary_table)