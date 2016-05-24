cv.MNS <-
function(dat, l1range, alpharange, K=5, parallel=FALSE, cores=NULL, verbose=FALSE){
  # wrapper function for CV 
  #
  # INPUT:
  #      - Dat: data used to fit MNS model
  #      - l1range: vector of potential l1 parameters
  #      - alpharange: vector of alpha parameters 
  #      - K: number of CV folds to perform
  #      - parallel, cores: should CV be run in parallel. If true then units indicates number of units/cores/workers
  #      - verbose: indicates whether progress should be printed. Only if parallel=FALSE
  #
  # OUTPUT:
  #      - l1: selected l1 value
  #      - alpha: selected alpha value
  #      - CV: grid of CV errors
  #
  
  if (parallel){
    return(CVparallel(Dat = dat, l1range = l1range, alpharange = alpharange, K = K, nodes=cores))
  } else {
    return(CVsequential(Dat = dat, l1range = l1range, alpharange = alpharange, K = K, verbose=verbose))
  }
  
  
}
