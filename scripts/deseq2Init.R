library("DESeq2")
library("readr")

# ---- Error Handling ----
if (!file.exists(snakemake@params[["sample_table"]])) {
  stop("Sample table not found: ", snakemake@params[["sample_table"]])
}

# ---- Read and Validate Sample Data ----
sample_table <- tryCatch(
  {
    tbl <- read_tsv(snakemake@params[["sample_table"]], 
                   col_types = cols(
                     sample = col_character(),
                     condition = col_factor(),
                     fq1 = col_character(),
                     fq2 = col_character()
                   ))
    
    # Validate required columns
    required_cols <- c("sample", "condition")
    if (!all(required_cols %in% colnames(tbl))) {
      stop("Sample table missing required columns: ", 
           paste(setdiff(required_cols, colnames(tbl)), collapse = ", "))
    }
    
    # Check for NA values
    if (any(is.na(tbl$condition))) {
      stop("NA values found in condition column")
    }
    
    tbl
  },
  error = function(e) {
    stop("Failed to read sample table: ", e$message)
  }
)

# ---- Verify Count Files Exist ----
count_files <- file.path("quantification/Deseq2Input", 
                        paste0(sample_table$sample, ".txt"))
missing_files <- count_files[!file.exists(count_files)]
if (length(missing_files) > 0) {
  stop("Missing count files:\n", paste(missing_files, collapse = "\n"))
}

# ---- Create DESeqDataSet ----
dds <- tryCatch(
  {
    # Create from HTSeqCount outputs
    dds <- DESeqDataSetFromHTSeqCount(
      sampleTable = sample_table,
      directory = "quantification/Deseq2Input",
      design = ~ condition
    )
    
    # Set reference level
    dds$condition <- relevel(dds$condition, ref = "untreated")
    dds
  },
  error = function(e) {
    stop("DESeqDataSet creation failed: ", e$message)
  }
)

# ---- Filter Low Count Genes ----
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# ---- Run DESeq2 ----
message("Running DESeq2 with ", ncol(dds), " samples and ", nrow(dds), " genes")
dds <- tryCatch(
  DESeq(dds, parallel = snakemake@threads > 1),
  error = function(e) {
    stop("DESeq2 analysis failed: ", e$message)
  }
)

# ---- Save Results ----
saveRDS(dds, file = snakemake@output[[1]])
message("DESeq2 object saved to ", snakemake@output[[1]])
