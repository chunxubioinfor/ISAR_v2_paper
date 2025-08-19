## Download the data which Jeroen suggested from https://zenodo.org/records/6826603 under Performance_GTEx.zip
## wget https://zenodo.org/records/6826603/files/Performance_GTEx.zip

# GTEx_data_prepare.R
# - Homogeneous GTEx subset
# - Generate 5v5 simulated DTU dataset (1 repeat)
# - NO pre-ISAR filtering (we'll filter inside ISAR later)
# - Save .Rdata with counts + design + truth

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(doMC)
})

# ---------------- Helper: NO-FILTER benchmark generator ----------------
getBenchmark_data_noFilter <- function(countData, metaData, nrRepList, fracGenesAffected) {
  lapply(seq_along(nrRepList), function(x) {
    
    # Step 1: Extract random sub-sample of correct size
    set.seed(x)
    localSampleSize <- nrRepList[[x]]
    
    if(TRUE) {
      localSubset <- sample(
        colnames(countData),
        localSampleSize * 2
      )
      
      localDesign <- data.frame(
        sample_id = localSubset,
        condition <- c(rep("a", localSampleSize), rep("b", localSampleSize)),
        stringsAsFactors = FALSE
      )
      localDesign <- localDesign[sort.list(localDesign$condition),]
    }
    
    # NO FILTERING: keep the selected samples as-is
    localCm <- as.matrix(countData[, localDesign$sample_id, drop = FALSE])
    
    # Build localTx from metaData (all transcripts present in localCm)
    localTx <- metaData[metaData$isoform_id %in% rownames(localCm), , drop = FALSE]
    
    # We still need to ensure genes chosen for swapping have >=2 isoforms
    gene_iso_counts <- table(localTx$gene_id)
    multi_genes <- names(gene_iso_counts)[gene_iso_counts >= 2]
    localTx <- localTx[localTx$gene_id %in% multi_genes, , drop = FALSE]
    localCm <- localCm[rownames(localCm) %in% localTx$isoform_id, , drop = FALSE]
    
    # Step 3: introduce DTU by swapping isoform abundances in condition B
    if(TRUE) {
      genesToModify <- sample(
        x = unique(localTx$gene_id),
        size = round(
          length(unique(localTx$gene_id)) * fracGenesAffected
        )
      )
      
      samplesToModify <- localDesign$sample_id[which(
        localDesign$condition == 'b')]
      
      cm_swapped <- localCm
      
      transcripts_toSwap_current <- c()
      transcripts_swapped_current <- c()
      transcripts_swapped_all <- c()
      transcripts_toSwap_all <- c()
      
      for (gene in genesToModify) {
        
        current <- localTx[which(localTx$gene_id==gene),]
        nSwap <- max(2,rbinom(1,nrow(current),1/3))
        
        transcripts_toSwap_current <- sample(
          x = current$isoform_id,
          size = nSwap)
        
        # swap order of txs completely by putting the fist one last
        transcripts_swapped_current <- c(transcripts_toSwap_current[-1],transcripts_toSwap_current[1])
        
        # Add to swapping queue
        transcripts_swapped_all <- c(transcripts_swapped_all,transcripts_swapped_current)
        transcripts_toSwap_all <- c(transcripts_toSwap_all,transcripts_toSwap_current)
      }
      
      # Perform the swapping in matrix
      cm_swapped[transcripts_toSwap_all,which(colnames(cm_swapped)%in%samplesToModify)] <- cm_swapped[transcripts_swapped_all,which(colnames(cm_swapped)%in%samplesToModify)]
      
      # Set swapping in localTx
      localTx$txSwapped <- vector(length=nrow(localTx))
      localTx[which(localTx$isoform_id %in% transcripts_swapped_all),"txSwapped"] <- TRUE 
      localTx$nrSamplesPerCondition <- localSampleSize
    }
    colnames(localTx) <- c('TXNAME','GENEID','txSwapped','nrSamplesPerCondition')
    
    # Combine data
    dataList <- list(
      data     = cm_swapped,
      design   = localDesign,
      metaInfo = localTx
    )
    
    return(dataList)
  })
}

# ---------------- Paths ----------------
data_dir      <- "~/Desktop/ISAR/ISAR_benchmark/data/GTEx"
counts_file   <- file.path(data_dir, "GTEx_counts.gz")     # renamed GTEx transcript expected counts
metadata_file <- file.path(data_dir, "GTEx_metadata.txt")  # renamed GTEx sample attributes

# ---------------- Load GTEx ----------------
message("Reading counts + metadata ...")
gtexCm     <- fread(counts_file, data.table = FALSE)
gtexSample <- fread(metadata_file, data.table = FALSE)

# ---------------- Subset to homogeneous samples ----------------
txInfo <- gtexCm[, c("transcript_id","gene_id")]
colnames(txInfo) <- c("isoform_id","gene_id")
rownames(gtexCm) <- txInfo$isoform_id
message("Subsetting samples to homogeneous set (Adrenal / Paxgene-plate / RNASEQ / B1) ...")
# Extract samples
gtexSample <- gtexSample[which(gtexSample$SMTSD == 'Adrenal Gland'),]
# Specifically filter on extraction kit
gtexSample <- gtexSample[which(gtexSample$SMNABTCHT == 'RNA Extraction from Paxgene-derived Lysate Plate Based'),]
# Specifically filter on SMAFRZE
gtexSample <- gtexSample[which(gtexSample$SMAFRZE == 'RNASEQ'),]
# Specifically filter on center
gtexSample <- gtexSample[which(gtexSample$SMCENTER == 'B1'),]
gtexCm <- gtexCm[,which(
  colnames(gtexCm) %in% gtexSample$SAMPID
)]

# ---------------- Simulation settings ----------------
samplesPrCondition <- 5L  # 5 vs 5
nrRepsMade         <- 1L  # one repeat
fracGenesAffected  <- 0.15
nrCoresToUse       <- 1L

if (nrCoresToUse > 1L) registerDoMC(cores = nrCoresToUse)

# Build the list structure expected by getBenchmark_data_* (single entry)
nrRepList <- setNames(
  list(rep(samplesPrCondition, nrRepsMade)),
  nm = paste0("samples_used_", samplesPrCondition, "_rep_1")
)

# ---------------- Generate dataset (no pre-filter) ----------------
message("Generating 5v5 dataset WITHOUT pre-filtering ...")
gtexBenchmarkData_NoFilter <- getBenchmark_data_noFilter(
  countData = gtexCm,
  metaData  = txInfo,
  nrRepList = nrRepList,
  fracGenesAffected = fracGenesAffected
)
names(gtexBenchmarkData_NoFilter) <- paste0(names(nrRepList), "_noFilter")


# ---------------- Save ----------------
out_rdata <- file.path(data_dir, "GTEx_benchmark_5v5_noFilter.Rdata")
save(gtexBenchmarkData_NoFilter, file = out_rdata)
message("Saved: ", out_rdata)


# ---------------- Quick summary ----------------
obj <- gtexBenchmarkData_NoFilter[[1]]
cat("\nSummary (first/only replicate):\n",
    "samples per condition: ", samplesPrCondition, "\n",
    "total samples: ", ncol(obj$data), "\n",
    "isoforms: ", nrow(obj$data), "\n",
    "multi-isoform truth genes affected (~15%): ",
    length(unique(obj$metaInfo$GENEID[obj$metaInfo$txSwapped])), " genes\n",
    sep = "")


