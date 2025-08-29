### Helper function for p-value calibration
pvcEstimateSigma2 <- function(
    ### The implementation of the p-value correction found in Ferguson et al P-value calibration for multiple testing problems in genomics ###
  # This is a modified version of the R function published with the article. The modification was simply to split the function into two functions.
  # This modified version only contains the parts for esimating and returning the Sigma2.
  # R code: http://www.math.bas.bg/~palejev/pvc/pvc.r
  # vignette : http://www.math.bas.bg/~palejev/pvc/
  
  ### Arguments
  pvalues,        # the original p-values
  pvaluesUse=NA,  # indicates which p-values would be used when estimating the likelihood. Possible values: either NA (default) or a boolean vector of the same length as pvalues. If NA, all p-values would be used when estimating the likelihood.If boolean vector of the same length as pvalues, only the elements of pvalues sliced by pvalueUse would be used when estimating the likelihood
  startsigma2=1,  # starting value of sigma2 in conditional likelihood and EM algorithm.
  condlik = TRUE, # set to FALSE if EM-algorithm based calibration is desired
  # conditional likelihood
  neff=NA,        # a prior estimate of the number of effective genes. Only useful when plotting the conditional likelihood
  plotl=FALSE,    # A logic indicting if a plot of the likelihood is desired, when running the conditional likelihood method
  le=0.1,          # lower bound for inclusion of p-values into the likelihood when running conditional likelihood method
  ue=1,           # upper bound for inclusion of p-values into the likelihood when running conditional likelihood method
  # EM implementation
  starteps=0.01,  # starting value for epsilon parameter in EM algorithm
  startmu=1,      # starting value for mu parameter in EM algorithm
  tol=.00001      # convergence criterion for EM algorithm
  
  ### Output
  # sigma2 - estimate of sigma2
) {
  Np <- length(pvalues)
  if((length(pvaluesUse)==1) && is.na(pvaluesUse)) {pvaluesUse <- rep(TRUE, Np)}
  if(length(pvaluesUse)!=Np) return("pvaluesUse not specified correctly")
  myindexes <- c(1:length(pvalues))[!is.na(pvalues)]
  pvalues <- pvalues[myindexes]
  pvaluesUse <- pvaluesUse[myindexes]
  prob_cal_needed <- NA
  
  ### Conditional likelyhood
  if(condlik){
    indexes <- c(1:length(pvalues))[pvalues < ue & pvalues >= le & pvaluesUse]
    tv <- qnorm(1-pvalues[indexes]/2)
    N <- length(tv)
    c1 <- qnorm(1-le/2)
    c2 <- qnorm(1-ue/2)
    simple_l <- function(sigma2){-(N/2)*log(sigma2)+sum(-1*tv^2/(2*sigma2))-N*log(pnorm(c1/sqrt(sigma2))-pnorm(c2/sqrt(sigma2)))}
    if(!is.na(neff)) simple_l <- function(sigma2){-(neff/2)*log(sigma2)+(neff/N)*sum(-1*tv^2/(2*sigma2))-neff*log(pnorm(c1/sqrt(sigma2))-pnorm(c2/sqrt(sigma2)))}
    sigma2new <- optim(par=startsigma2,method="L-BFGS-B",function(x){-1*simple_l(x)},lower=0.01, upper=10)$par
    if(plotl) {
      sigma_possible <- seq(from = 0.2, to = 2, by = 0.0001)
      like_vals <- sapply(sigma_possible, function(x){simple_l(x)})
      like_vals <- like_vals - max(like_vals)
      cum_area <- cumsum(exp(like_vals))/sum(exp(like_vals))
      prob_cal_needed <- max(cum_area[sigma_possible < 0.9])
      plot(sigma_possible, exp(like_vals)/(.0001*sum(exp(like_vals))),type="l",main="", xlab=expression(sigma^2),ylab = "")
      title(expression(paste("Normalized likelihood for ",sigma^2)))
      S <- .0001*sum(exp(like_vals))
      for(j in 1:length(cum_area[sigma_possible < 0.9])) lines(x=rep(sigma_possible[j],2),y=c(0,exp(like_vals[j])/S),col="palevioletred")
      text(x=.75,y=.5,labels=paste("area = ",round(prob_cal_needed,2),sep=""))
    }
  }
  
  ### EM algorithm
  if(!condlik){
    indexes <- c(1:length(pvalues))[pvalues<1 & pvaluesUse]
    tv <- qnorm(1-pvalues[indexes]/2)
    tv[tv=="Inf"]=max(tv[tv!="Inf"])
    N <- length(tv)
    
    sigma2old <- startsigma2
    currenterror <- 1
    muold <- startmu
    epsold <- starteps
    while(currenterror > tol){
      ddiff <- dnorm((tv-muold)/sqrt(sigma2old))+dnorm((-1*tv-muold)/sqrt(sigma2old))
      dnormal <- 2*dnorm(tv/sqrt(sigma2old))
      probs <- epsold*ddiff/((epsold*ddiff)+(1-epsold)*dnormal)
      a1 <- epsold*dnorm((tv-muold)/sqrt(sigma2old))/((epsold*ddiff)+(1-epsold)*dnormal)
      a2 <- epsold*dnorm((-1*tv-muold)/sqrt(sigma2old))/((epsold*ddiff)+(1-epsold)*dnormal)
      munew <- weighted.mean(c(tv,-1*tv), w=c(a1,a2))
      sigma2new <- (1/N)*sum((1-probs)*tv^2 + a1*(tv-munew)^2+a2*(-1*tv-munew)^2)
      epsnew <- mean(probs)
      currenterror <- max(abs(epsnew-epsold),abs(sigma2new-sigma2old),abs(munew-muold))
      sigma2old <- sigma2new
      muold <- munew
      epsold <- epsnew
      
    }
    postprobs= rep(NA,length(pvalues))
    postprobs[indexes]=probs
    if(length(myindexes)<Np){
      postprobs1 <- rep(NA,Np)
      postprobs1[myindexes] <- postprobs
      postprobs <- postprobs1
    }
    logL <- -1*sum((1-probs)*tv^2)/(2*sigma2new)-1*sum(a1*(tv-munew)^2)/(2*sigma2new)-1*sum(a2*(-1*tv-munew)^2)/(2*sigma2new)-(N/2)*log(2*pi*sigma2new)+log(epsnew)*sum(probs)+log(1-epsnew)*(N-sum(probs))
  }
  
  return(sigma2new)
  
}

pvcApplySigma <- function(
    ### The implementation of the p-value correction found in Ferguson et al P-value calibration for multiple testing problems in genomics ###
  # This is a modified version of the R function published with the article. The modification was simply to split the function into two functions.
  # This modified version only contains the parts for using the estimated sigma2 to callibrated the p-values.
  # R code: http://www.math.bas.bg/~palejev/pvc/pvc.r
  # vignette : http://www.math.bas.bg/~palejev/pvc/
  
  ### Arguments
  pvalues,        # the original p-values
  sigma2new
  
  ### Output
  # calibrated p-values
) {
  Np <- length(pvalues)
  myindexes <- c(1:length(pvalues))[!is.na(pvalues)]
  pvalues <- pvalues[myindexes]
  
  ### transform the original pvalues.
  indexes <- c(1:length(pvalues))[pvalues<1]
  tv <- qnorm(1-pvalues[indexes]/2)
  tv[tv=="Inf"]=max(tv[tv!="Inf"])
  N <- length(tv)
  pvalues= rep(1,length(pvalues))
  pvalues[indexes] <- pchisq(tv^2/sigma2new,0,lower.tail=FALSE,df=1)
  if(length(myindexes)<Np){
    pvalues1 <- rep(NA,Np)
    pvalues1[myindexes] <- pvalues
    pvalues <- pvalues1
  }
  
  return(pvalues)
}

# map ENST -> ENSG from obj
truth_map <- obj$metaInfo %>% select(TXNAME, GENEID) %>% distinct()

# isoformFeatures minimal table (ENST, ENSG)
feat <- bm_isar_v2_saturn$isoformFeatures %>% as.data.frame() %>%
  select(isoform_id, gene_id)
# overwrite gene_id with ENSG if available
feat <- feat %>%
  left_join(truth_map, by = c("isoform_id"="TXNAME")) %>%
  mutate(gene_id = ifelse(!is.na(GENEID), GENEID, gene_id)) %>%
  select(-GENEID)


# start from the object
expr_df <- bm_isar_v2_saturn$isoformRepExpression %>% as.data.frame()
# move 'isoform_id' into rownames and drop the column
stopifnot("isoform_id" %in% colnames(expr_df))
rownames(expr_df) <- expr_df$isoform_id
expr_df$isoform_id <- NULL
# now a clean matrix: rows = ENST, cols = samples
expr_mat <- as.matrix(expr_df)

design <- bm_isar_v2_saturn$designMatrix %>% as.data.frame() %>%
  transmute(sampleID = sampleID, group = as.character(condition))
grpA <- design$sampleID[design$group=="a"]; nA <- length(grpA)
grpB <- design$sampleID[design$group=="b"]; nB <- length(grpB)

mean_se <- function(v){
  n <- sum(!is.na(v))
  c(mean = mean(v, na.rm = TRUE),
    se   = if (n > 1) sd(v, na.rm = TRUE)/sqrt(n) else NA_real_)
}
# isoform-level
iso_tbl <- tibble(isoform_id = rownames(expr_mat)) %>%
  mutate(
    iso_value_1  = apply(expr_mat[, grpA, drop=FALSE], 1, function(x) mean_se(x)["mean"]),
    iso_stderr_1 = apply(expr_mat[, grpA, drop=FALSE], 1, function(x) mean_se(x)["se"]),
    iso_value_2  = apply(expr_mat[, grpB, drop=FALSE], 1, function(x) mean_se(x)["mean"]),
    iso_stderr_2 = apply(expr_mat[, grpB, drop=FALSE], 1, function(x) mean_se(x)["se"])
  ) %>%
  left_join(feat, by="isoform_id")

# gene-level per-sample sums, then mean/SE
sum_by_gene <- function(cols) {
  as.data.frame(expr_mat[, cols, drop=FALSE]) %>%
    tibble::rownames_to_column("isoform_id") %>%
    left_join(feat, by="isoform_id") %>%
    select(-isoform_id) %>%
    group_by(gene_id) %>%
    summarize(across(everything(), sum, na.rm=TRUE), .groups="drop") %>%
    tibble::column_to_rownames("gene_id")
}
gene_expr_A <- sum_by_gene(grpA)
gene_expr_B <- sum_by_gene(grpB)
gene_ids <- union(rownames(gene_expr_A), rownames(gene_expr_B))

gene_tbl <- tibble(gene_id = gene_ids) %>%
  mutate(
    gene_value_1  = ifelse(gene_id %in% rownames(gene_expr_A),
                           apply(gene_expr_A[gene_id[gene_id %in% rownames(gene_expr_A)], , drop=FALSE], 1, function(x) mean_se(x)["mean"]),
                           NA_real_),
    gene_stderr_1 = ifelse(gene_id %in% rownames(gene_expr_A),
                           apply(gene_expr_A[gene_id[gene_id %in% rownames(gene_expr_A)], , drop=FALSE], 1, function(x) mean_se(x)["se"]),
                           NA_real_),
    gene_value_2  = ifelse(gene_id %in% rownames(gene_expr_B),
                           apply(gene_expr_B[gene_id[gene_id %in% rownames(gene_expr_B)], , drop=FALSE], 1, function(x) mean_se(x)["mean"]),
                           NA_real_),
    gene_stderr_2 = ifelse(gene_id %in% rownames(gene_expr_B),
                           apply(gene_expr_B[gene_id[gene_id %in% rownames(gene_expr_B)], , drop=FALSE], 1, function(x) mean_se(x)["se"]),
                           NA_real_)
  )

# merge isoform + gene stats, compute IF and variances
df <- iso_tbl %>%
  left_join(gene_tbl, by="gene_id") %>%
  mutate(
    IF1 = ifelse(gene_value_1 > 0, iso_value_1 / gene_value_1, NA_real_),
    IF2 = ifelse(gene_value_2 > 0, iso_value_2 / gene_value_2, NA_real_),
    dIF = IF2 - IF1,
    nrReplicates_1 = nA,
    nrReplicates_2 = nB,
    gene_var_1 = (gene_stderr_1 * sqrt(nrReplicates_1))^2,
    gene_var_2 = (gene_stderr_2 * sqrt(nrReplicates_2))^2,
    iso_var_1  = (iso_stderr_1  * sqrt(nrReplicates_1))^2,
    iso_var_2  = (iso_stderr_2  * sqrt(nrReplicates_2))^2,
    IF_var_1 = (1 / (gene_value_1^2)) * ( iso_var_1 + (IF1^2 * gene_var_1) - (2 * IF1 * iso_var_1) ),
    IF_var_2 = (1 / (gene_value_2^2)) * ( iso_var_2 + (IF2^2 * gene_var_2) - (2 * IF2 * iso_var_2) )
  )

# -------- Welch t-test on IF difference --------
welch_df <- function(v1, n1, v2, n2) {
  num <- (v1/n1 + v2/n2)^2
  den <- ifelse(n1>1, (v1^2)/(n1^2*(n1-1)), NA_real_) + ifelse(n2>1,(v2^2)/(n2^2*(n2-1)), NA_real_)
  num/den
}

res <- df %>%
  mutate(
    valid = is.finite(IF1) & is.finite(IF2) &
      is.finite(IF_var_1) & is.finite(IF_var_2) &
      IF_var_1 > 0 & IF_var_2 > 0,
    meanDiff = IF2 - IF1,
    seDiff   = sqrt( (IF_var_1 / nrReplicates_1) + (IF_var_2 / nrReplicates_2) ),
    tscore   = ifelse(valid, meanDiff / seDiff, NA_real_),
    df_welch = ifelse(valid, welch_df(IF_var_1, nrReplicates_1, IF_var_2, nrReplicates_2), NA_real_),
    pValue   = ifelse(valid, 2*pt(abs(tscore), df=df_welch, lower.tail=FALSE), NA_real_)
  )

# -------- Calibrate and BH adjust --------
# Use highly expressed isoforms for sigma^2 estimation (as in v1)
isoExpCutoff <- max(1, stats::quantile(c(res$iso_value_1, res$iso_value_2), probs=0.50, na.rm=TRUE))
he <- res %>% filter(iso_value_1 > isoExpCutoff, iso_value_2 > isoExpCutoff, !is.na(pValue))
sigma2 <- pvcEstimateSigma2(he$pValue)  # condlik=TRUE by default

calibrated <- if (!is.na(sigma2) && sigma2 < 0.9) pvcApplySigma(res$pValue, sigma2) else res$pValue
q_iso <- p.adjust(calibrated, method = "BH")

res2 <- res %>%
  mutate(isoform_switch_q_value = q_iso)

# -------- Gene-level FDR: min isoform q per gene --------
gene_table <- res2 %>%
  group_by(gene_id) %>%
  summarise(gene_switch_q_value = suppressWarnings(min(isoform_switch_q_value, na.rm = TRUE)),
            .groups = "drop")

x_grid <- seq(0.005, 0.25, by = 0.005)
curve_v1 <- tibble(
  thresh = x_grid,
  n_sig  = map_int(x_grid, ~ sum(gene_table$gene_switch_q_value <= .x, na.rm = TRUE)),
  version = "ISAR v1 (manual)"
)
intended_fdr <- 0.05
p <- ggplot(curve_v1, aes(x = thresh, y = n_sig)) +
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

gene_table <- gene_table %>% 
  mutate(is_truth = gene_id %in% truth_genes)

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