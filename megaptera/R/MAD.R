MAD <- function(DNAbin, model = "JC69"){
  
  if ( nrow(DNAbin) > 1 ){
    p.uncorr <- dist.dna(DNAbin, model = "raw", pairwise.deletion = TRUE)
    p.corr <- dist.dna(DNAbin, model = model, pairwise.deletion = TRUE) 
    x <- mad(p.corr - p.uncorr, na.rm = TRUE)
  } else {
    x <- 0
  }
  x
}