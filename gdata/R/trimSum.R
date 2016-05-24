### trimSum.R
###------------------------------------------------------------------------
### What: Sum trimmed values - code
### $Id$
### Time-stamp: <2008-12-20 12:11:27 ggorjan>
###------------------------------------------------------------------------

trimSum <- function(x, n, right=TRUE, na.rm=FALSE, ...)
{
  ## --- Setup ---

  if(!is.vector(x) | is.list(x))
    stop("'x' must be a vector - for now")
  if(!is.numeric(x))
    stop("'x' must be numeric")
  if(length(x) <= n)
    stop("'n' must be smaller than the length of x")
  
  ## --- Trim ---

  N <- length(x)
  if(right) {
    x2 <- x[1:n]
    x2[n] <- sum(x[n:N], na.rm=na.rm)
  } else {
    k <- (N - n + 1)
    x2 <- x[k:N]
    x2[1] <- sum(x[1:k], na.rm=na.rm)
  }
  
  ## --- Return ---

  x2
}

###------------------------------------------------------------------------
### trimSum.R ends here
