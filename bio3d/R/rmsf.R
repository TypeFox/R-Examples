"rmsf" <-
function(xyz, average = FALSE) {
  if(is.null(dim(xyz)))
    stop("input 'xyz' has NULL dimension")

  ## Cov function changed ~ R.2.7
  my.sd <- function (x, na.rm = FALSE) {
    if (is.matrix(x))
      apply(x, 2, my.sd, na.rm = na.rm)
    else if (is.vector(x)) {
      if(na.rm) x <- x[!is.na(x)]
      if(length(x) == 0) return(NA)
      sqrt(var(x, na.rm = na.rm))
    }
    else if (is.data.frame(x))
      sapply(x, my.sd, na.rm = na.rm)
    else {
      x <- as.vector(x)
      my.sd(x,na.rm=na.rm)
    }
  }

  fluct = rowSums( matrix( my.sd(xyz, na.rm = TRUE), ncol=3, byrow=TRUE )^2, na.rm=TRUE )
  if(average) {
     if(ncol(xyz) %% 3 == 0)
        d = ncol(xyz) / 3
     else 
        d = ncol(xyz)
     return( sqrt( sum(fluct)/d ) )
  } else {
     return( sqrt( fluct ) )
  }
  
}

