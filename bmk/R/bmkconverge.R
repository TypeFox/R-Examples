#'  MCMC chain convergence diagnostic.
#'  
#'  This takes an MCMC chain and divides it into batches of size \code{binsize} and calculates the Hellinger distance between consecutive batches.
#'  
#'  @title bmkcoverge:  Convergence via the Hellinger distance
#'  @param inputlist1 A list of the MCMC chains
#'  @param binsize a scalar giving how large each bin should be for consecutive batches.
#'  outputs the Hellinger distances between the sampled distribution for one scenario against the other.
#'  @references Boone EL, Merrick JR and Krachey MJ.  
#'   A Hellinger distance approach to MCMC diagnostics.
#'   Journal of Statistical Computation and Simulation, 
#'   \code{DOI:10.1080/00949655.2012.729588}.
#'  @export
#'  @examples
#'  \dontrun{
#' library(dismo); library(MCMCpack); 
#' data(Anguilla_train)
#' b0mean <- 0
#' b0precision <- (1/5)^2
#' mcmclen = 1000
#' burn=10000
#' MCMC.one <- MCMClogit(Angaus ~ SegSumT+DSDist+USNative+as.factor(Method)+DSMaxSlope+USSlope, 
#'                 data=Anguilla_train,burnin=burn, mcmc=mcmclen, beta.start=-1, 
#'                 b0=b0mean, B0=b0precision)
#'  }
#'  data(MCMCsamples)
#'  mcmclen <- 1000
#'  bmkconverge(MCMC.one,mcmclen/10)
#' 
bmkconverge = function(inputlist1, binsize=1000){
  require(functional)
  require(plyr)
  require(coda)
  
  uncolnames = Compose(colnames, unique)
  lenuni = Compose(unique, length)
  
  n.variables = ncol(inputlist1)
  n.iter = nrow(inputlist1)
  n.var = ncol(inputlist1)
  ########################################################################
  
  dataset1 = inputlist1
  size = (floor(n.iter/binsize)-1)
  output = matrix(0, nrow = n.var, ncol =size)
  for(i in 1:n.variables){
    print(length(HWconverg1(dataset1[,i], binsize)))
    output[i,] <- HWconverg1(dataset1[,i], binsize)
  }
  rownames(output) <- colnames(dataset1)
  colnames(output) <- (1:size)*binsize
  output
}
