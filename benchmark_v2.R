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

# ------------------ Paths ------------------
data_dir <- "~/Desktop/ISAR/ISAR_benchmark/data/GTEx"
out_dir  <- "~/Desktop/ISAR/ISAR_benchmark/v2_out"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
gtf_file <- file.path(data_dir, "gencode.v26.chr_patch_hapl_scaff.annotation.gtf.gz")

# ------------------ Load dataset ------------------
load(file.path(data_dir, "GTEx_benchmark_5v5_noFilter.Rdata"))  # loads gtexBenchmarkData_NoFilter
obj <- gtexBenchmarkData_NoFilter[[1]]

# ------------------ Build inputs ------------------
isoformCountMatrix <- obj$data %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "isoform_id") %>%     # REQUIRED first column
  dplyr::relocate(isoform_id)
designMatrix <- obj$design %>%
  dplyr::transmute(
    sampleID  = sample_id,                               # REQUIRED name
    condition = as.character(condition)                  # strings
  )

# -------- Build SwitchAnalyzeRlist --------
bm_isar_v2 <- importRdata(
  isoformCountMatrix   = isoformCountMatrix,
  designMatrix         = designMatrix,
  isoformExonAnno      = gtf_file
)

# Set all the arguments to 0, just filter the single-isoform genes to avoid crash of saturn
bm_isar_v2_filtered <- preFilter(
  switchAnalyzeRlist = bm_isar_v2,
  isoCount = 0,
  min.Count.prop = 0,
  IFcutoff = 0,
  min.IF.prop = 0,
  alpha = 0.05,
  dIFcutoff = 0.1,
  quiet = FALSE
)


bm_isar_v2_saturn <- isoformSwitchTestSatuRn(
  bm_isar_v2_filtered,
  alpha = 0.05,
  dIFcutoff = 0.1,
  reduceToSwitchingGenes = FALSE
)

extractSwitchSummary(bm_isar_v2_saturn)

#> extractSwitchSummary(bm_isar_v2_saturn)
#Comparison nrIsoforms nrSwitches nrGenes
#1     a vs b       2443       2344    1898


# Extract gene-level FDR directly
gene_table <- as_tibble(bm_isar_v2_saturn$isoformFeatures) %>%
  select(gene_id, gene_switch_q_value) %>%
  distinct()   # one row per gene

# Build #sig vs FDR curve
x_grid <- seq(0.005, 0.25, by = 0.005)
curve_v2 <- tibble(
  thresh = x_grid,
  n_sig  = map_int(x_grid, ~ sum(gene_table$gene_switch_q_value <= .x, na.rm = TRUE)),
  version = "ISAR v2"
)

# Plot
intended_fdr <- 0.05
p <- ggplot(curve_v2, aes(x = thresh, y = n_sig)) +
  geom_point() +
  geom_vline(xintercept = intended_fdr, linetype = "dashed") +
  labs(
    title    = "Comparison of current and previous IsoformSwitchAnalyzeR",
    subtitle = "Number of significant genes as a function of False Discovery Rate.\nThe intended FDR is highlighted by the dashed line.",
    x = "False Discovery Rate threshold",
    y = "Number of significant genes"
  ) +
  theme_bw(base_size = 12)
p

ggsave(file.path(out_dir, "scatter_v2_genelevel.png"), p, width = 6.5, height = 4.8, dpi = 300)

# TPR vs FDR curve
truth_map <- obj$metaInfo %>%
  select(TXNAME, GENEID) %>%
  distinct()
bm_isar_v2_saturn$isoformFeatures <- bm_isar_v2_saturn$isoformFeatures %>%
  as_tibble() %>%
  left_join(truth_map, by = c("isoform_id" = "TXNAME")) %>%
  mutate(gene_id = ifelse(!is.na(GENEID), GENEID, gene_id)) %>%
  select(-GENEID) %>%
  as.data.frame() 

truth_genes <- unique(obj$metaInfo$GENEID[obj$metaInfo$txSwapped == TRUE])

gene_table <- as_tibble(bm_isar_v2_saturn$isoformFeatures) %>%
  select(gene_id, gene_switch_q_value) %>%
  distinct() %>%
  mutate(is_truth = gene_id %in% truth_genes)

total_truth <- sum(gene_table$is_truth)

x_grid <- seq(0.005, 0.25, by = 0.005)

# total truth genes
total_truth <- sum(gene_table$is_truth)

x_grid <- seq(0.005, 0.25, by = 0.005)

curve_tpr <- tibble(thresh = x_grid) %>%
  mutate(
    TP = map_int(thresh, ~ {
      sum(gene_table$is_truth[!is.na(gene_table$gene_switch_q_value) &
                                gene_table$gene_switch_q_value <= .x])
    }),
    FP = map_int(thresh, ~ {
      sum(!gene_table$is_truth[!is.na(gene_table$gene_switch_q_value) &
                                 gene_table$gene_switch_q_value <= .x])
    })
  ) %>%
  mutate(
    TPR    = TP / total_truth,
    empFDR = ifelse((TP+FP) == 0, NA_real_, FP / (TP+FP))
  )

# Plot TPR vs thresholded FDR
ggplot(curve_tpr, aes(x = thresh, y = TPR)) +
  geom_line() +
  geom_vline(xintercept = 0.05, linetype = "dashed") +
  labs(
    title = "TPR vs FDR",
    x = "Nominal FDR threshold",
    y = "True Positive Rate"
  ) +
  theme_bw()


