## This function provides a matrix with all combinations of a vector K
## (containing degrees) and vector I (0/1 for exclude/include for
## factors).

matrix.combn <- function(K.vec1,K.vec2=NULL, num.x=0,num.z=NULL) {

  if(num.x==0 && num.z==0) stop(" must provide at least one variable")

  ls <- list()
  if(num.x>0) {
      for(i in 1:num.x) ls[[i]] <- K.vec1
      if(!is.null(K.vec2)) {
          for(i in 1:num.x) ls[[num.x+i]] <- K.vec2
          if(!is.null(num.z)) for(i in 1:num.z) ls[[2*num.x+i]] <- 0:1
      } else {
          if(!is.null(num.z)) for(i in 1:num.z) ls[[num.x+i]] <- 0:1
      }
  }
  else 
  {
      if(!is.null(num.z)) for(i in 1:num.z) ls[[num.x+i]] <- 0:1
  }
  return(as.matrix(do.call(expand.grid,ls)))

}
