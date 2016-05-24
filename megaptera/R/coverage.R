# CALCULATE COVERAGE of species and nucleotides in alignments
# package: megaptera
# called by: stepD, filter.alignment
# author: Christoph Heibl
# last change 2014-10-29

coverage <- function(DNAbin, what = "species"){
  
  if ( !inherits(DNAbin, "DNAbin") )
    stop("'DNAbin' must be of class 'DNAbin'")
  if ( !is.matrix(DNAbin) )
    stop("'DNAbin' must be alignend")
  
  what <- match.arg(what, c("species", "nucleotides"))
  
  b <- as.raw(c(136, 40, 72, 24))
  percentInformation <- function(x, b){
    length(which(x %in% b))
  }
  if ( what == "species" ){
    m <- apply(DNAbin, 1, percentInformation, b)
    m <- m/ncol(DNAbin)
  } else {
    m <- apply(DNAbin, 2, percentInformation, b)
    m <- m/nrow(DNAbin)
  }
  m
}