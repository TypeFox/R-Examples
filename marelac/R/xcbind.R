## non-exported function called in diffcoeff.R
## wrapper for cbind with additional checks of arguments
##
## extensions:
## - checks if one or both arguments are NULL or have length zero
## - converts list or vector to data frame
## - works if length of both matches after recyling rule
## - otherwise cbind throws error
## limitation: accepts only two arguments
## author: thpe

xcbind <- function(x, y) {
  makedf <- function(x) {
    if (is.list(x)) {
      x <- as.data.frame(x)
    } else if (is.vector(x)) {
      x <- data.frame(x)
    }
    x
  }
  
  ## length=0 if x is NULL or numeric(0)
  if (length(x) & length(y)) return(cbind(makedf(x), makedf(y)))
  if (length(y)) return(makedf(y))
  if (length(x)) return(makedf(x))
  return(NULL)
}
