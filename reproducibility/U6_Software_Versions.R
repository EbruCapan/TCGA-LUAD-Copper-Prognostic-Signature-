# U6 - Record current software environment; does not prove historical versions.
packages <- c('TCGAbiolinks', 'DESeq2', 'edgeR', 'limma', 'survival',
              'survminer', 'glmnet', 'timeROC', 'ConsensusClusterPlus',
              'ComplexHeatmap', 'maftools', 'clusterProfiler',
              'oncoPredict', 'openxlsx', 'ggplot2')
outdir <- file.path('results', 'U6_Reproducibility')
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
versions <- data.frame(
  Package = packages,
  Installed = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
  Version = vapply(packages, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p))
    else NA_character_
  }, character(1))
)
write.csv(versions, file.path(outdir, 'U6_Software_Package_Versions.csv'),
          row.names = FALSE)
capture.output(sessionInfo(), file = file.path(outdir, 'U6_R_SessionInfo.txt'))
print(versions, row.names = FALSE)
cat('R: ', R.version.string, '\n')
