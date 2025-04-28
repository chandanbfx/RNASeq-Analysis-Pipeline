library(ggplot2)
library(svglite)
library(fontquiver)

# Input validation
if (!file.exists(snakemake@input[[1]])) {
  stop("Input results file not found: ", snakemake@input[[1]])
}

# Read data with proper headers
res <- tryCatch(
  read.table(
    file = snakemake@input[[1]], 
    header = TRUE, 
    row.names = 1,
    sep = "\t",
    check.names = FALSE
  ),
  error = function(e) stop("Failed to read results file: ", e$message)
)

# Validate required columns
required_cols <- c("log2FoldChange", "pvalue", "padj")
if (!all(required_cols %in% colnames(res))) {
  stop("Missing required columns in results. Needed: ", 
       paste(required_cols, collapse = ", "))
}

# Font setup
fonts <- fontquiver::font_families("Liberation")
fonts$symbol$symbol <- fontquiver::font_symbol("Symbola")

# Create plot
svglite(snakemake@output[["vol_plot"]], user_fonts = fonts)

tryCatch({
  alpha <- 0.05
  lfc_threshold <- 1.5  # Configurable threshold
  
  # Enhanced coloring
  res$significant <- ifelse(
    abs(res$log2FoldChange) > lfc_threshold & res$padj < alpha,
    "Significant",
    "Not significant"
  )
  
  # Create base plot
  plot(
    res$log2FoldChange, 
    -log10(res$padj),
    col = ifelse(res$significant == "Significant", "red", "gray50"),
    pch = 20,
    cex = 0.8,
    panel.first = grid(),
    main = "Volcano Plot", 
    xlab = expression(log[2]("Fold Change")),
    ylab = expression(-log[10]("Adjusted p-value"))
  
  # Add threshold lines
  abline(v = 0, lty = 2)
  abline(v = c(-lfc_threshold, lfc_threshold), col = "blue", lty = 2)
  abline(h = -log10(alpha), col = "blue", lty = 2)
  
  # Label top significant genes
  sig_genes <- res[res$significant == "Significant", ]
  if (nrow(sig_genes) > 0) {
    top_genes <- head(sig_genes[order(sig_genes$padj), ], 10)  # Label top 10
    text(
      top_genes$log2FoldChange,
      -log10(top_genes$padj),
      labels = rownames(top_genes),
      pos = 3,
      cex = 0.6,
      col = "black"
    )
  }
  
  # Add legend
  legend(
    "topright",
    legend = c("Significant", "Not significant"),
    col = c("red", "gray50"),
    pch = 20,
    bty = "n"
  )
}, error = function(e) {
  dev.off()
  stop("Volcano plot generation failed: ", e$message)
})

dev.off()
