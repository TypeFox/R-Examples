# Function to calculate Importance Sampling (IS) weights from (log) kernel and 
# (log) candidate density evaluations at M parameter values
#  
# inputs:
#    lnkernel : [vector size M] of log- kernel density evaluations 
#    lncand   : [vector size M] of log- candidate density evaluations
# outputs:
#    r        : [vector size M] of IS weights
#
# author : Nalan Basturk
# date   : 20120912
fn.ISwgts <- function(lnkernel,lncand){
  if(!is.vector(lnkernel) | !is.vector(lncand))
    stop("inputs of 'fn.ISwgts' should be vectors")
  if(length(lnkernel)!=length(lncand))
    stop("inputs of 'fn.ISwgts' should be vectors of the same size")
  r <- lnkernel - lncand
  r <- r-max(r)   # robustification
  r <- exp(r)
  r  = r/sum(r)   # weights sum up to 1
  return(r)
}