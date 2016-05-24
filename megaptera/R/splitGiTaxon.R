splitGiTaxon <- function(x){
  x <- unlist(strsplit(x, "_"))
  c(gi = tail(x, 1), taxon = paste(head(x, -1), collapse = "_"))
}