daeEnv <- new.env()
tol <- c(.Machine$double.eps ^ 0.5, .Machine$double.eps ^ 0.5)
names(tol) <- c("element.tol", "eigen.tol")
assign("daeTolerance",  tol, envir=daeEnv)

#set.daeTolerance <- function(tolerance = daeTolerance)
#{ proj <- environment(set.daeTolerance)
#  unlockBinding("daeTolerance", proj)
#  daeTolerance <<- tolerance
#  lockBinding("daeTolerance", proj)
#  invisible(daeTolerance)
#}

set.daeTolerance <- function(element.tol=NULL, eigen.tol=NULL)
{ tol <- get("daeTolerance", envir=daeEnv)
  
  #Set value for element.tol, if one supplied
  if (!is.null(element.tol))
  { if (!is.numeric(element.tol))
      stop("non-numeric argument")
    if (length(element.tol) != 1)
    { element.tol <- element.tol[1]
      warning("More than one element.tol value supplied - only first used")
    } 
    tol[1] <- element.tol
  }
    
  #Set value for eiggen.tol, if one supplied
  if (!is.null(eigen.tol))
  { if (!is.numeric(eigen.tol))
      stop("non-numeric argument")
    if (length(eigen.tol) != 1)
    { eigen.tol <- eigen.tol[1]
      warning("More than one eigen.tol value supplied - only first used")
    } 
    tol[2] <- eigen.tol
  }
  
  names(tol) <- c("element.tol", "eigen.tol")
  assign("daeTolerance", tol, envir = daeEnv)
  invisible(tol)
}

get.daeTolerance <- function()
{ tol <- get("daeTolerance", envir=daeEnv)
  return(tol)
}

#function to test whether all elements are zero
#was all(abs(x) < daeTolerance)
"is.allzero" <- function(x)
{ daeTolerance <- get("daeTolerance", envir=daeEnv)
  mean(abs(x), na.rm=TRUE) < daeTolerance[["element.tol"]]
}


remove.repeats <- function(x, tolerance = 1E-06)
  # function to remove repeated values that differ by no more than daeTolerance  
{ #daeTolerance <- get("daeTolerance", envir=daeEnv)
  n <- length(x)
  if (n > 1)
  { repeats <-   c(FALSE, abs(x[2:n] - x[1:(n-1)]) < tolerance)
    x <- x[!repeats]
  }
  return(x)
}

"check.arg.values" <- function(arg.val, options)
  #Function to check that arg.val is one of the allowed values
  #and to return the position of the argument in the set of values
  #that is stored in options
{ kopt <- pmatch(arg.val, options)
  if (is.na(kopt))
    stop("Value for argument, ",arg.val,", is not an allowed option")
  return(kopt)
}

"ginv" <- function(x, tol = .Machine$double.eps ^ 0.5)
{ 
  # computes Moore-Penrose inverse of a matrix
  if (!is.matrix(x) | length(dim(x)) != 2 )
    stop("x must be a matrix")
  svd.x <- svd(x)
  nonzero.x <- (svd.x$d > svd.x$d[1] * tol)
  rank.x <- sum(nonzero.x)
  geninv.x <- matrix(0, dim(x)[1], dim(x)[2])
  if (rank.x)
  { i <- matrix((1:length(nonzero.x))[nonzero.x], rank.x, 2)
    geninv.x[i] <- 1/svd.x$d[nonzero.x]
    if (all(nonzero.x))
      geninv.x <- svd.x$v %*% geninv.x %*% t(svd.x$u)
    else 
      geninv.x <- svd.x$v[, nonzero.x] %*% geninv.x[nonzero.x, nonzero.x] %*% 
      t(svd.x$u[, nonzero.x])
  }
  geninv.x
}
