#'  Hellinger distance between distributions
#'  
#'  This computes the Hellinger distance for all pairwise combinations of MCMC chains.
#'
#'  @param chains1 A matrix of MCMC for the same variable.  Each column corresponds to a different chain.
#'  @return c2 A vector containing the minimum and maximum Hellinger distances across all pairwise comparisons.
#'  @note The matrix must consist of samples from the same variable derived from different chains.
HBconverg1 <- function(chains1){
  chains2 <- as.matrix(chains1)
  n1 <- nrow(chains2)
  n2 <- ncol(chains2)
  c1 <- 0
  for(i in 1:(n2-1)){
    for(j in (i+1):(n2)){
      out1 <- HDistNoSize(chains2[,i],chains2[,j])
      c1 <- c(c1,out1)
    }
  }
  c2 <- cbind(min(c1[2:length(c1)]),max(c1[2:length(c1)]))
  return(c2)
}
