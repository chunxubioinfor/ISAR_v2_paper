truth_iso <- obj$metaInfo %>%
  transmute(
    isoform_id = TXNAME,
    gene_id    = GENEID,
    truth_swapped = txSwapped
  )
alpha <- 0.05
dIFcutoff <- 0.10

## 1) Isoform-level: every significant switching isoform
sig_iso <- extractTopSwitches(
  switchAnalyzeRlist = bm_isar_v2_analyzed_dexseq,
  filterForConsequences = FALSE, # set TRUE only if you've run analyzeSwitchConsequences()
  alpha = alpha,
  dIFcutoff = dIFcutoff,
  n = Inf                        # <- all of them
)
# columns include: gene_id, isoform_id, dIF, isoform_switch_p_value, isoform_switch_q_value, etc.

## 2) Gene-level: unique genes that have ≥1 significant switching isoform
sig_gen <- sig_iso %>%
  dplyr::distinct(gene_id) %>%
  dplyr::arrange(gene_id)

res <- bm_isar_v2_analyzed_dexseq$isoformFeatures
res <- bm_isar_v2_analyzed_saturn$isoformFeatures
## isoform-level
sig_iso_manual <- res %>%
  dplyr::filter(!is.na(isoform_switch_q_value),
                isoform_switch_q_value <= alpha,
                is.na(dIF) | abs(dIF) >= dIFcutoff)

## gene-level (gene is significant if any isoform is significant)
sig_gen_manual <- sig_iso_manual %>%
  dplyr::distinct(gene_id)

truth_map <- obj$metaInfo %>%
  dplyr::select(TXNAME, GENEID) %>%   # TXNAME = ENST, GENEID = ENSG
  dplyr::distinct()
isar_map <- bm_isar_v2_analyzed_dexseq$isoformFeatures %>%
  dplyr::select(isoform_id, gene_id) %>%
  dplyr::distinct()
isar_map <- bm_isar_v2_analyzed_saturn$isoformFeatures %>%
  dplyr::select(isoform_id, gene_id) %>%
  dplyr::distinct()
map_full <- isar_map %>%
  dplyr::left_join(truth_map, by = c("isoform_id" = "TXNAME"))
sig_iso_ENSG <- sig_iso_manual %>%
  dplyr::left_join(map_full, by = c("isoform_id" = "isoform_id"))
truth_genes <- unique(obj$metaInfo$GENEID[obj$metaInfo$txSwapped == TRUE])
sig_genes <- unique(sig_iso_ENSG$GENEID)
intersect(truth_genes,sig_genes)

