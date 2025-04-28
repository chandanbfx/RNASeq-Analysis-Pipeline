library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("svglite")
library("fontquiver")

# Validate input
if (!file.exists(snakemake@input[[1]])) {
  stop("Input RDS file not found: ", snakemake@input[[1]])
}

# Load data
dds <- tryCatch(
  readRDS(snakemake@input[[1]]),
  error = function(e) stop("Failed to read DESeq2 object: ", e$message)
)

# Font setup
fonts <- fontquiver::font_families("Liberation")
fonts$symbol$symbol <- fontquiver::font_symbol("Symbola")

# Create heatmap
svglite(snakemake@output[[1]], user_fonts = fonts)

tryCatch({
  vsd <- vst(dds, blind = FALSE)
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  
  # Improved labeling
  rownames(sampleDistMatrix) <- paste(vsd$condition, 
                                     colnames(vsd), 
                                     sep = " - ")
  colnames(sampleDistMatrix) <- NULL
  
  # Color palette
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  pheatmap(
    sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col = colors,
    fontsize_row = 8,
    fontsize_col = 8
  )
}, error = function(e) {
  dev.off()
  stop("Heatmap generation failed: ", e$message)
})

dev.off()
