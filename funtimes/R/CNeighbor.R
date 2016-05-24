CNeighbor <- function(Bu, Bv, Alpha, Beta, Delta, Theta){
  p <- length(Bu)
  if(is.null(dim(Bv)[2])) {Bv <- matrix(Bv, ncol=1)}
  colSums(abs(Bu-Bv)/(Beta-Alpha) <= Delta) >= Theta*p
}