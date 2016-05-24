### getXY.R

getXY <- function(x, y = NULL, unidim.allowed = TRUE)
  ## Author: Rene Locher
  ## Version: 2006-06-20
{
  if (missing(x)) xarg <- "x" else xarg <- deparse(substitute(x))

  if (missing(y)) yarg <- "y" else yarg <- deparse(substitute(y))

  if (is.matrix(x) | is.data.frame(x)) {
    if (ncol(x)>1) {
      if (is.null(y)) y <- x[,2] else
      stop("'",xarg,
           "' must have only 1 column when y is supplied separately\n")
    }
    x <- x[,1]
  } else if (is.list(x)) {
    if (length(x)>1) {
      y <- x[[2]]
      x <- x[[1]]
      if (length(y)!=length(x))
        stop("First and second element of the list must have identical length!\n")
    } else x <- x[[1]]
  }

  if (is.null(y)) {
    if (!unidim.allowed) stop("'",yarg,"' is not defined!\n")
    y <- x
    x <- 1:length(x)
  } else {
    y <-  unlist(y)
    if (length(y)!=length(x))
      stop("Vector '",yarg,"' and 'x' must have identical lengths!\n")
  }

 return(data.frame(x = x, y = y))
} ## getXY

