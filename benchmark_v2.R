# ------------------ Download ------------------
if (!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools")
}
if (!requireNamespace("pfamAnalyzeR", quietly = TRUE)){
  devtools::install_github("kvittingseerup/pfamAnalyzeR")
}

if (!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools")
}
devtools::install_github("kvittingseerup/IsoformSwitchAnalyzeR", build_vignettes = TRUE)

suppressPackageStartupMessages({
  library(IsoformSwitchAnalyzeR)   # v2
  library(tidyverse)
})

# ------------------ Paths (edit) ------------------
data_dir <- "~/Desktop/ISAR/ISAR_benchmark/data/GTEx"
out_dir  <- "~/Desktop/ISAR/ISAR_benchmark/v2_out"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
gtf_file <- file.path(data_dir, "gencode.v26.chr_patch_hapl_scaff.annotation.gtf.gz")

# ------------------ Load dataset ------------------
load(file.path(data_dir, "GTEx_benchmark_5v5_noFilter.Rdata"))  # loads gtexBenchmarkData_NoFilter
obj <- gtexBenchmarkData_NoFilter[[1]]

# ------------------ Build inputs ------------------
isoformCountMatrix <- as.data.frame(obj$data)

designMatrix <- obj$design %>% as.data.frame()
rownames(designMatrix) <- designMatrix$sample_id
designMatrix$condition <- factor(designMatrix$condition, levels = c("a","b"))

# -------- Build SwitchAnalyzeRlist --------
aSL <- importRdata(
  isoformCountMatrix   = isoformCountMatrix,
  designMatrix         = designMatrix,
  isoformExonAnno      = gtf_file,
  removeNonConvensionalChr = FALSE,
  quiet = TRUE
)

