as.yearmon2 <- function(x, ...){
  clx <- class(x)
  {
    if((length(clx)==1) && (clx=="numeric") && (round(x, 2) == x)){
      ix <- format(x, nsmall=2)
      ox <- as.yearmon(ix, "%Y.%m")
      names(ox) <- ix 
    }
    else{
      ox <- as.yearmon(x, ...)
      names(ox) <- x
    }
  }
  if((nox <- length(unique(ox))) != (nx <- length(x))){
    warning(nx-nox, " duplicate months found in 'x'; ",
            " returning 'x' unchanged")
    return(x)
  }
  ox 
}

