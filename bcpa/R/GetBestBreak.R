#' Find most likely change point in irregular time series
#'
#' Finds the single best change point according to the likelihood function.  Used internally within \code{\link{WindowSweep}}.
#'
#' @param x  vector of time series values.
#' @param t	vector of times of measurements associated with x.
#' @param range tange of possible breaks. Default (0.6) runs approximately from 1/5th to 4/5ths of the total length of the time series.
#' @param ... additional parameters to pass to \code{\link{GetDoubleL}} function.

#' @return returns a single row (vector) with elements: \code{breaks},\code{tbreaks},\code{mu1},\code{sigma1},\code{rho1},\code{LL1},\code{mu2},\code{sigma2},\code{rho2},\code{LL2},\code{LL}. The breakpoint is calculated for a range of possible values of width \code{range*l} (where \code{l} is the length of the time series). The output of this function feeds \code{\link{WindowSweep}}.
#' 
#' @seealso  \code{\link{WindowSweep}} which uses it, and \code{\link{GetDoubleL}} for the likelihood estimation. 
#' @author Eliezer Gurarie
#' 
#' @examples 
#' # An example with a single break:
#' x <- c(arima.sim(list(ar = 0.9), 20) + 10, arima.sim(list(ar = 0.1), 20)) 
#' t <- 1:length(x)
#' plot(t,x, type="l")
#' (bb <- GetBestBreak(x,t, tau=FALSE))
#' abline(v = bb[2], col=2)

GetBestBreak <-
  function(x,t,range=0.6, ...)
  {
    lower <- round((1-range)/2 * length(t))
    
    GetDoubleL2 <- function(breaks, ...)
      # get the total likelihood
      sum(GetDoubleL(x,t,breaks, ...)[c(4,8)])
    
    BestBreak <- round(optimize(GetDoubleL2, lower = lower, upper=length(t) - lower, maximum=TRUE)$max)
    
    myDoubleL <- GetDoubleL(x,t,BestBreak, ...)
    c(bb.index = BestBreak, bb.time = t[BestBreak], myDoubleL, ll.total = myDoubleL[4] + myDoubleL[8])
  }