deleteGaps <- function(x, nmax = nrow(x) - 4){
  
  if (!inherits(x, "DNAbin")) stop("'x' is not of class 'DNAbin'")
  
	gapCounter <- function(n, nmax){
	  length(which(n %in% as.raw(4))) > nmax
	}
  id <- apply(x, 2, gapCounter, nmax = nmax)
  x[, !id]
}