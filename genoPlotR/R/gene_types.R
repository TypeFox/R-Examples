################################################################################
# Plot helpers: list and determine gene "types"
################################################################################
gene_types <- function(auto=TRUE){
  types <- c("arrows", "headless_arrows",
             "blocks", "bars", "points", "text", "lines", 
             "side_blocks", "side_bars", "side_points", "side_text",
             "side_lines",
             "introns", "exons", "side_exons")
  if (auto) types <- c("auto", types)
  types
}
auto_gene_type <- function(n_genes){
  if (max(n_genes) > 1000){
    gene_type <- "side_bars"
  } else if (max(n_genes) > 100) {
    gene_type <- "side_blocks"
  } else {
    gene_type <- "side_bars"
  }
  gene_type
}
