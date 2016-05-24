#' Plot distribution of CDF
#' 
#' Plot the elicited pointwise median and credible interval for an uncertain population CDF
#'  
#' @param medianfit The output of a \code{fitdist} command following elicitation
#'  of the expert's beliefs about the population median.
#' @param precisionfit The output of a \code{fitdist} command following elicitation
#'  of the expert's beliefs about the population precision.
#' @param lower lower limit on the x-axis for plotting.
#' @param upper upper limit on the x-axis for plotting.
#' @param ql lower quantile for the plotted pointwise credible interval. 
#' @param qu upper quantile for the plotted pointwise credible interval.
#' @param median.dist The fitted distribution for the population median. Can be one of \code{"normal"},
#'  \code{"lognormal"} or \code{"best"}, where \code{"best"} will select the best fitting out of 
#'  normal and lognormal.
#' @param precision.dist The fitted distribution for the population precision. Can either be \code{"gamma"}
#'  or \code{"lognormal"}. 
#' @param n.rep The number of randomly sampled CDFs used to estimated the median
#'  and credible interval.
#' @param n.X The number of points on the x-axis at which the CDF is evaluated.
#' @param fontsize Font size used in the plots.
#'  
#' @examples \dontrun{
#' prfit <- fitprecision(interval = c(60, 70), propvals = c(0.2, 0.4), trans = "log")
#' medianfit <- fitdist(vals = c(50, 60, 70), probs = c(0.05, 0.5,  0.95), lower = 0)
#' cdfplot(medianfit, prfit)
#'  }
#' @import ggplot2
#' @export

cdfplot <- function(medianfit, precisionfit, 
                    lower = NA, upper = NA, ql = 0.025, 
                    qu = 0.975, median.dist = "best", precision.dist = "gamma",
                    n.rep = 10000, n.X = 100, fontsize = 18){
  
  if(precision.dist!="gamma" & precision.dist!="lognormal"){
    stop('precision.dist must equal one of "gamma" or "lognormal"')
  }
  
  ps <- psample(medianfit, precisionfit, lower, 
                upper, median.dist, precision.dist, n.rep, n.X)
  X <- ps$X
  pX <- ps$pX
  pX <- apply(pX, 2, sort)
  
  medcdf <- lowercdf <- uppercdf <- NULL # hack to avoid R CMD check NOTE
  
  plotdata <- data.frame(X=X, medcdf=pX[ceiling(0.5 * n.rep),],
                         lowercdf=pX[ceiling(ql * n.rep), ], 
                         uppercdf=pX[ceiling(qu * n.rep), ])
  theme_set(theme_grey(base_size = fontsize))
  cdp <- ggplot(plotdata, aes(x = X, y = medcdf)) +
    ylab("P(X<=x)") +
    xlab("x") +
    geom_ribbon(aes(ymin = lowercdf, 
                    ymax = uppercdf), colour = "blue", fill = "blue", 
                alpha=0.5) +
    geom_line(size = 2.0, colour = "red") 
  print(cdp)
}