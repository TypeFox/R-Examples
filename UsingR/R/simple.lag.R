##' simple lag plot
##'
##' From code posted by Martyn Plummer \code{<plummer "at" iarc.fr>}
##' @param x data
##' @param lag value of lab
##' @param FUN function to apply
##' @return NULL
##'
##' @export
"simple.lag" <-
  function(x, lag, FUN=mean)
{
  ## from Martyn Plummer <plummer@iarc.fr>
  FUN=match.fun(FUN)
  n <- length(x)
  y <- numeric(n-lag)
  for (i in (1+lag):n) {
    y[i - lag] <- FUN(x[(i-lag):i])
  }
  return(y)
}
