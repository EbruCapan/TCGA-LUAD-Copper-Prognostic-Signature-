# U5 - Time-dependent ROC for the six-gene TCGA-LUAD prognostic risk score
# Run from the project root (LUAD_Copper_Project).
# Input is a patient-level RDS generated upstream; this script does NOT rebuild the risk score.

required <- c('timeROC', 'openxlsx')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop('Install required packages: ', paste(missing, collapse = ', '))

input <- file.path('drug_sensitivity', 'results', 'U2_6Gene_504',
                   'TCGA_504_Final_Risk_Groups.rds')
if (!file.exists(input)) stop('Input not found: ', input, '\nRun from project root.')
risk504 <- readRDS(input)
needed <- c('submitter_id', 'OS.time', 'OS', 'risk_score')
stopifnot(all(needed %in% names(risk504)), nrow(risk504) == 504L,
          !anyDuplicated(risk504$submitter_id),
          !anyNA(risk504[, needed]), all(risk504$OS %in% c(0, 1)),
          all(is.finite(risk504$OS.time)), all(risk504$OS.time > 0),
          all(is.finite(risk504$risk_score)))

outdir <- file.path('results', 'U5_Time_Dependent_ROC_504')
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
years <- c(1, 3, 5)
times <- years * 365.25
roc <- timeROC::timeROC(T = risk504$OS.time, delta = risk504$OS,
                        marker = risk504$risk_score, cause = 1,
                        weighting = 'marginal', times = times, iid = TRUE)
stopifnot(length(roc$AUC) == length(years))
ci <- stats::confint(roc, level = 0.95)$CI_AUC / 100
results <- data.frame(
  Year = years, Time_Days = times, AUC = as.numeric(roc$AUC),
  CI_Lower_95 = as.numeric(ci[, 1]), CI_Upper_95 = as.numeric(ci[, 2]),
  Events_By_Time = vapply(times, function(t)
    sum(risk504$OS == 1 & risk504$OS.time <= t), integer(1)),
  At_Risk = vapply(times, function(t)
    sum(risk504$OS.time >= t), integer(1))
)
print(results)
saveRDS(roc, file.path(outdir, 'U5_TimeROC_6Gene_504.rds'))
saveRDS(results, file.path(outdir, 'U5_AUC_Results_504.rds'))

colors <- c('#2166AC', '#E69F00', '#B2182B')
draw_roc <- function() {
  par(mar = c(5, 5, 2.5, 2), family = 'sans')
  plot(roc, time = times[1], col = colors[1], lwd = 2.8, title = FALSE)
  plot(roc, time = times[2], add = TRUE, col = colors[2], lwd = 2.8)
  plot(roc, time = times[3], add = TRUE, col = colors[3], lwd = 2.8)
  abline(a = 0, b = 1, col = 'gray65', lty = 2, lwd = 1.4)
  title(main = 'Time-dependent ROC analysis', cex.main = 1.2)
  labels <- sprintf('%d-year AUC = %.3f (95%% CI %.3f-%.3f)',
                    results$Year, results$AUC,
                    results$CI_Lower_95, results$CI_Upper_95)
  legend('bottomright', legend = labels, col = colors, lwd = 2.8,
         bty = 'n', cex = 0.85)
}
pdf(file.path(outdir, 'U5_Figure6_TimeROC_FINAL.pdf'),
    width = 8.5, height = 7.5, useDingbats = FALSE)
draw_roc()
dev.off()
png(file.path(outdir, 'U5_Figure6_TimeROC_FINAL_600dpi.png'),
    width = 8.5, height = 7.5, units = 'in', res = 600, type = 'cairo')
draw_roc()
dev.off()

info <- data.frame(
  Parameter = c('Cohort', 'Patients', 'Deaths', 'Censored',
                'Prognostic genes', 'Outcome', 'Method', 'Validation'),
  Value = c('TCGA-LUAD', nrow(risk504), sum(risk504$OS == 1),
            sum(risk504$OS == 0),
            'RHOV, MT1H, LOXL2, SOD3, SLC31A2, MT1A',
            'Overall survival', 'Time-dependent ROC, marginal IPCW',
            'Apparent performance in development cohort')
)
file <- file.path(outdir, 'U5_TimeROC_Results_Reproduced.xlsx')
openxlsx::write.xlsx(list(Time_Dependent_AUC = results,
                          Analysis_Information = info),
                     file = file, overwrite = TRUE)
check <- openxlsx::read.xlsx(file, sheet = 'Time_Dependent_AUC')
stopifnot(nrow(check) == 3L, all(c('AUC', 'CI_Lower_95', 'CI_Upper_95') %in% names(check)))
cat('U5 complete. Outputs: ', normalizePath(outdir), '\n')
