library("DESeq2")
library("svglite")
library("fontquiver")

# Error handling for input
if (!file.exists(snakemake@input[[1]])) {
  stop("Input RDS file not found: ", snakemake@input[[1]])
}

# Load data
dds <- tryCatch(
  readRDS(snakemake@input[[1]]),
  error = function(e) stop("Failed to read RDS file: ", e$message)
)

# Font setup
fonts <- fontquiver::font_families("Liberation")
fonts$symbol$symbol <- fontquiver::font_symbol("Symbola")

# Create plot
svglite(snakemake@output[[1]], user_fonts = fonts)

tryCatch({
  counts <- rlog(dds, blind = FALSE)
  plotPCA(counts, intgroup = snakemake@params[["pca_labels"]])
}, error = function(e) {
  dev.off()
  stop("PCA plotting failed: ", e$message)
})

dev.off()
