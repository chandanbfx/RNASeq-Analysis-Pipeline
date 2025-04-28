library("DESeq2")
library("apeglm")
library("ggplot2")

# ---- Input Validation ----
if (!file.exists(snakemake@input[[1]])) {
  stop("Input RDS file not found: ", snakemake@input[[1]])
}

# ---- Load Data ----
dds <- tryCatch(
  {
    dds <- readRDS(snakemake@input[[1]])
    # Validate object type
    if (!inherits(dds, "DESeqDataSet")) {
      stop("Input is not a DESeqDataSet object")
    }
    dds
  },
  error = function(e) {
    stop("Failed to load DESeq2 object: ", e$message)
  }
)

# ---- Validate Contrast ----
contrast <- snakemake@params[["contrast"]]
if (length(contrast) != 2) {
  stop("Contrast must contain exactly 2 conditions")
}
if (!all(contrast %in% levels(dds$condition))) {
  stop("Invalid contrast conditions. Available levels: ", 
       paste(levels(dds$condition), collapse = ", "))
}

# ---- Run Differential Expression ----
message("\nRunning DE analysis for contrast: ", paste(contrast, collapse = " vs "))

res <- tryCatch(
  {
    # Initial results
    res <- results(
      dds,
      contrast = c("condition", contrast[1], contrast[2]),
      alpha = 0.05,
      independentFiltering = TRUE
    )
    
    # LFC shrinkage
    coef_name <- resultsNames(dds)[grep(paste0("condition_", contrast[1], "_vs_", contrast[2]), 
                                      resultsNames(dds))]
    if (length(coef_name) == 1) {
      resLFC <- lfcShrink(
        dds,
        coef = coef_name,
        type = "apeglm",
        res = res
      )
      # Transfer shrunk values
      res$log2FoldChange <- resLFC$log2FoldChange
      res$lfcSE <- resLFC$lfcSE
    } else {
      warning("Could not find matching coefficient for LFC shrinkage")
    }
    
    # Sort by adjusted p-value
    res[order(res$padj), ]
  },
  error = function(e) {
    stop("DE analysis failed: ", e$message)
  }
)

# ---- MA Plot ----
tryCatch(
  {
    svg(snakemake@output[["ma_plot"]])
    plotMA(res, 
           main = paste("MA Plot:", paste(contrast, collapse = " vs ")),
           ylim = c(-3, 3),
           alpha = 0.05)
    dev.off()
  },
  error = function(e) {
    if (exists("svg_device")) dev.off()
    stop("MA plot generation failed: ", e$message)
  }
)

# ---- Save Results ----
write.table(
  as.data.frame(res),
  file = snakemake@output[["table"]],
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

message("\nResults saved to:\n",
        "- Table: ", snakemake@output[["table"]], "\n",
        "- MA Plot: ", snakemake@output[["ma_plot"]])
