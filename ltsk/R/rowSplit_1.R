rowSplit <- function(x,f,tran=F)
{
  ## x: a matrix to be divided into groups
  ## f: a 'factor' defines the groups
  ## tran: default split by row
  if(tran) {
    stopifnot(ncol(x)==length(f)) 
  }
  else{
    stopifnot(nrow(x)==length(f)) 
  }
  if(tran) x <- t(x)
  idx <- split(1:nrow(x),f)
  #browser()
  if(tran){
    r <- lapply(idx,function(elmt) t(x[elmt,,drop=FALSE]))
  }
  else{
    r <- lapply(idx,function(elmt) x[elmt,,drop=FALSE])    
  }
  r
}