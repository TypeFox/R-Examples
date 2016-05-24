#' running function (default median absolute deviation)
#' @param aref data array
#' @param k window in data points, default 300
#' @param func default med but can be any function taking a vector and returning a summary
#' @examples
#' x = rnorm(500)
#' x = c(x,rnorm(1000,3,2))
#' x = c(x,runif(1000,4,6))
#' y = runFun(x,k=51,func=mad)
#' hist(y)#[500:490]
#' y2 = runFun(x,k=51,func=median)
#' plot(x,pch="*")
#' lines(y2,col=2,lwd=3)
#' lines(y2+y,col=3,lwd=3)
#' lines(y2-y,col=3,lwd=3)
#' tic = runFun(x,k=51,func=function(x,...){mean(x)})
#' plot(x,pch=".")
#' abline(h=0,col=2)
#' lines(tic,col=3,lwd=3)
#' @export
#' @seealso  \code{\link{runmed}}
runFun <- function(aref,k=301,func=mad)
{
  m = k %/% 2
  N = length(aref)
  res = rep(NA,N)
  for(i in 1:(N-k)){
    sub = aref[i:(i+k)]
    res[i + m]=func(sub,na.rm=TRUE)
  }
  res[1:m] = rep(res[m+1],m)
  res[(N-m ):N] = rep(res[N-m-1],m+1)
  return(res)
}
