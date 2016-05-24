# generate null sp x site matrices with same total abundance as observed and species probabilities drawn abundances within the base null matrix
generate_nul_resample <- function(nul, obs)
{
  nul_resample <- matrix(0, nrow=nrow(nul),ncol=ncol(nul),dimnames=dimnames(nul)) # auto-set sparseness depending on filling
  for (i in 1:ncol(nul_resample))
  {
    # fill in resamples with the same total abundance
    samples <- sample(1:nrow(nul),size=sum(obs[,i]),replace=TRUE,prob=nul[,i])
    nul_resample[,i] <- tabulate(samples, nbins=nrow(nul))    
  }	
  
  return(nul_resample)
}  

