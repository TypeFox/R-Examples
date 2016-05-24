# ========================================================================
# cbind  -  column bind for rvs
# ========================================================================
# Note: It's inconvenient that we cannot call the generic (default) cbind
#       if class attributes are set.
#

# DEBUG: cbind(1, rvnorm(1), rvnorm(1)) causes an error msg
#   In cbind(v, unclass(x[[i]]), deparse.level = deparse.level) :
#     number of rows of result is not a multiple of vector length (arg 2)

cbind.rv <- function(..., deparse.level = 1)
{
  if (deparse.level != 1) 
    .NotYetUsed("deparse.level != 1")
  x <- list(...)
  if (length(x)<1) return(NULL)
  v <- NULL
  for (i in seq(along=x)) {
    v <- cbind(v, unclass(x[[i]]), deparse.level=deparse.level)
  }
  class(v) <- class(rv())
  return(v)
}


## ========================================================================
## rvbind.rv  -  row bind for rvs
## ========================================================================
##

rbind.rv <- function(..., deparse.level = 1)
{
  if (deparse.level != 1) 
    .NotYetUsed("deparse.level != 1")
  x <- list(...)
  if (length(x)<1) return(NULL)
  v <- NULL
  for (i in seq(along=x)) {
    v <- rbind(v, unclass(x[[i]]), deparse.level=deparse.level)
  }
  class(v) <- class(rv())
  return(v)
}

