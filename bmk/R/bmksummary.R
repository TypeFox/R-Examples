#'  Summaries of MCMC chains
#'  
#'  Generate the mean, standard deviation, median, 2.5 percent and 97.5 percent quantiles, Gelman Rubin statistic for convergence, effective samples size and the minimum and maximum Hellinger distances across all chains.  outputs summaries for the MCMC samples including the convergence diagnostics of Gelman and Rubin and the Hellinger distance of Boone, Merrick and Krachey.
#'  
#'
#'  @param inputlist A list of the combined MCMC chains for all samples from one scenario.
#'  @references Boone EL, Merrick JR and Krachey MJ.  
#'   A Hellinger distance approach to MCMC diagnostics.
#'   Journal of Statistical Computation and Simulation, 
#'   \code{DOI:10.1080/00949655.2012.729588}.
#'  @export
#'  @examples
#'  data(MCMCsamples)
#'  bmksummary(list( MCMC.one, MCMC.two, MCMC.three ))
##  \dontrun{
##  library(dismo); library(MCMCpack)
##  data(Anguilla_train)
##  b0mean <- 0
##  b0precision <- (1/5)^2
##  mcmclen = 1000
##  burn=10000
##  MCMC.one <- MCMClogit(Angaus ~~ SegSumT+DSDist+USNative+as.factor(Method)+DSMaxSlope+USSlope, 
##                   data=Anguilla_train,burnin=burn, mcmc=mcmclen, beta.start=-1, 
##                   b0=b0mean, B0=b0precision)
##  MCMC.two <- MCMClogit(Angaus ~~ SegSumT+DSDist+USNative+as.factor(Method)+DSMaxSlope+USSlope, 
##                   data=Anguilla_train,burnin=burn, mcmc=mcmclen, beta.start=-.5, 
##                   b0=b0mean, B0=b0precision)
##  MCMC.three <- MCMClogit(Angaus ~~ SegSumT+DSDist+USNative+as.factor(Method)+DSMaxSlope+USSlope, 
##                     data=Anguilla_train,burnin=burn, mcmc=mcmclen, beta.start=0, 
##                     b0=b0mean, B0=b0precision)
##  bmksummary( list( MCMC.one, MCMC.two, MCMC.three ) ) 
##  }
bmksummary = function(inputlist){
  require(functional)
  require(plyr)
  require(coda)
  
  uncolnames = Compose(colnames, unique)
  lenuni = Compose(unique, length)
  
  
  n.variables = ncol(inputlist[[1]])
  n.chains = length(inputlist)
  n.iter = nrow(inputlist[[1]])
  ####
  CI = Curry(quantile, probs=c(0.025,0.975))
  catf = Curry(cat,fill=TRUE)
  pastes = Curry(paste, sep="")
  std.error = function(x) sqrt(var(x))/sqrt(length(x))
  n.var = ncol(inputlist[[1]])
  
  
  colfun= function(x,fun) fun(matrix(x, ncol=1))
  dataset = Reduce(rbind, inputlist)
  for(i in 1:n.variables){
    output = matrix(0,nrow = n.var, ncol = 10)
    for(j in 1:n.var){ 
      mcmc.obj = list(NA)
      colfun2 = Curry(colfun, x=dataset[,j])
      data = matrix(dataset[,j], nrow=n.iter, ncol=n.chains,byrow=FALSE)
      for(kk in 1:n.chains)mcmc.obj[[kk]] = mcmc(data[,kk])
      mcmc.obj = mcmc.list(mcmc.obj)
      output[j,c(1:3, 9:10)] = c(mean(dataset[,j]),median(dataset[,j]),std.error(dataset[,j]),CI(dataset[,j]))
      output[j,4] = effectiveSize(mcmc.obj)[[1]]
      output[j,5] <- n.iter
      output[j,6] = gelman.diag(mcmc.obj)$psrf[[1]]
      output[j,7:8] = HBconverg1(data)
    }
  }
  colnames(output) = c("Mean","Median","SE", "EffSS","Samples", "GRB", "BMKMin", "BMKMax" ,"0.025CI","0.975CI")
  rownames(output) = colnames(inputlist[[1]])
  output
}
