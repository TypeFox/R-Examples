#'  Hellinger distance within consecutive batches of MCMC samples.
#'  
#'  Determine if a specific chain has converged.  This takes a chain and divides it into batches and calculates the Hellinger distance between consecutive batches.
#'
#'  @param chain1 A vector of a single MCMC chain.
#'  @param batchsize1 An integer that defines the size of the batches.
#'  @return c2 A vector of Hellinger distances between consecutive batches.
HWconverg1 <- function(chain1, batchsize1 = 1000){
  chain2 <- as.matrix(chain1)
  batchlist1 <- seq(1,nrow(chain2),by=batchsize1)
  if(is.integer(length(chain1)/batchsize1)==FALSE){print("Warning: First Batch Will Be Shorter Than Rest Since Chain Length Is Not a Multiple of Batch Size")}
  c1 <- 0
  for(i in 1:(length(batchlist1)-1)){
    batchlabel1 <- 1:(i*batchsize1)
    batchlabel2 <- batchlist1[i+1]:((i+1)*batchsize1)
    out1 <- HDistNoSize(chain1[batchlabel1],chain1[batchlabel2])
    c1 <- c(c1,out1)
  }
  c2 <- c1[2:length(c1)]
  return(c2)
}
