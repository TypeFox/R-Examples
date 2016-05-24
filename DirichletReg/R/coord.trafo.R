coord.trafo <- function(Y){

  .Deprecated(msg="coord.trafo() is now deprecated.\nTo plot in the new (qua)ternary coordinate systems, use toSimplex()")
  
  return(toSimplex(Y[,c(2L, 3L, 1L)]))

  if((ncol(Y) < 3) | (ncol(Y) > 4)) stop("does not work for less than 3 or more than 4 dimensions")

  if(ncol(Y) == 3){
    x <- 0.5 * (Y[,3] + 2*Y[,2]) / rowSums(Y)
    y <- sqrt(3)/2 * Y[,3] / rowSums(Y)
    return(cbind(x,y))
  }

  if(ncol(Y) == 4){
    x <- 0.5 * (Y[,3] + 2*Y[,2]) / rowSums(Y) + 0.5 * Y[,4]
    y <- sqrt(3)/2 * Y[,3] / rowSums(Y) + sqrt(3)/6 * Y[,4] / rowSums(Y)
    z <- sqrt(3)/2 * Y[,4] / rowSums(Y)
    return(cbind(x,y,z))
  }

}
