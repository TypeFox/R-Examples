#'  Hellinger distance between two MCMC chains for sensitivity studies
#'
#'  Determine if two identically dimensioned sets of chains match.  This is good for conducting sensitivity studies.
#'  
#'  @param inputlist1 A list of the combined MCMC chains for all samples from one scenario.
#'  @param inputlist2 A list of the combined MCMC chains for all samples from another scenario.
#'  @references Boone EL, Merrick JR and Krachey MJ.  
#'   A Hellinger distance approach to MCMC diagnostics.
#'   Journal of Statistical Computation and Simulation, 
#'   \code{DOI:10.1080/00949655.2012.729588}.
#'  @export
#'  @examples
#'  data(MCMCsamples)
#'  bmksensitive(MCMC.one.mean0, MCMC.one.mean1)
#'  \dontrun{
#'  library(dismo); library(MCMCpack)
#'  data(Anguilla_train)
#'  b0mean0 <- 0
#'  b0mean1 <- 1
#'  b0precision <- (1/5)^2
#'  mcmclen = 1000
#'  burn=10000
#'  MCMC.one.mean0 <- MCMClogit(Angaus ~ SegSumT+DSDist+USNative+as.factor(Method)+DSMaxSlope+USSlope, 
#'                   data=Anguilla_train,burnin=burn, mcmc=mcmclen, beta.start=-1, 
#'                   b0=b0mean0, B0=b0precision)
#'  MCMC.one.mean1 <- MCMClogit(Angaus ~ SegSumT+DSDist+USNative+as.factor(Method)+DSMaxSlope+USSlope, 
#'                   data=Anguilla_train,burnin=burn, mcmc=mcmclen, beta.start=-.5, 
#'                   b0=b0mean1, B0=b0precision)
#'  bmksensitive(one, two)
#'  }
bmksensitive = function(inputlist1, inputlist2){
  require(functional)
  require(plyr)
  
  uncolnames = Compose(colnames, unique)
  lenuni = Compose(unique, length)
  
  n.variables = ncol(inputlist1)
  n.iter = nrow(inputlist1)
  n.var = ncol(inputlist1)
  dataset1 = inputlist1
  dataset2 = inputlist2
  output = matrix(0, nrow = n.var, ncol =1 )
  for(i in 1:n.variables){
    output[i,1] <- HDistNoSize(dataset1[,i],dataset2[,i])
  }
  
  rownames(output) <- colnames(dataset1)
  colnames(output) <- "HellingerDist"
  output
}
