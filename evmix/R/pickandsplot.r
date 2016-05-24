#' @export
#' 
#' @title Pickands Plot
#'
#' @description Produces the Pickand's plot.
#'
#' @inheritParams hillplot
#' 
#' @details Produces the Pickand's plot including confidence intervals.
#'   
#'   For an ordered iid sequence \eqn{X_{(1)}\ge X_{(2)}\ge\cdots\ge X_{(n)}} 
#'   the Pickand's estimator of the reciprocal of the shape parameter \eqn{\xi} 
#'   at the \eqn{k}th order statistic is given by 
#'   \deqn{\hat{\xi}_{k,n}=\frac{1}{\log(2)} \log\left(\frac{X_{(k)}-X_{(2k)}}{X_{(2k)}-X_{(4k)}}\right).}
#'   Unlike the Hill estimator it does not assume positive data, is valid for any \eqn{\xi} and
#'   is location and scale invariant.
#'   The Pickands estimator is defined on orders \eqn{k=1, \ldots, \lfloor n/4\rfloor}. 
#'   
#'   Once a sufficiently low order statistic is reached the Pickand's estimator will
#'   be constant, upto sample uncertainty, for regularly varying tails. Pickand's
#'   plot is a plot of \deqn{\hat{\xi}_{k,n}} against the \eqn{k}. Symmetric asymptotic
#'   normal confidence intervals assuming Pareto tails are provided.
#'   
#'   The Pickand's estimator is for the GPD shape \eqn{\xi}, or the reciprocal of the
#'   tail index \eqn{\alpha=1/\xi}. The shape is plotted by default using
#'   \code{y.alpha=FALSE} and the tail index is plotted when \code{y.alpha=TRUE}.
#'   
#'   A pre-chosen threshold (or more than one) can be given in
#'   \code{try.thresh}. The estimated parameter (\eqn{\xi} or \eqn{\alpha}) at
#'   each threshold are plot by a horizontal solid line for all higher thresholds. 
#'   The threshold should be set as low as possible, so a dashed line is shown
#'   below the pre-chosen threshold. If Pickand's estimator is similar to the
#'   dashed line then a lower threshold may be chosen.
#'   
#'   If no order statistic (or threshold) limits are provided 
#'   \code{orderlim = tlim = NULL} then the lowest order statistic is set to \eqn{X_{(1)}} and
#'   highest possible value \eqn{X_{\lfloor n/4\rfloor}}. However, Pickand's estimator is always
#'   output for all \eqn{k=1, \ldots, \lfloor n/4\rfloor}.
#'   
#'   The missing (\code{NA} and \code{NaN}) and non-finite values are ignored.
#'   
#'   The lower x-axis is the order \eqn{k}. The upper axis is for the corresponding threshold.
#' 
#' @return \code{\link[evmix:pickandsplot]{pickandsplot}} gives Pickand's plot. It also 
#'   returns a dataframe containing columns of the order statistics, order, Pickand's
#'   estimator, it's standard devation and \eqn{100(1 - \alpha)\%} confidence
#'   interval (when requested).
#' 
#' @note 
#' Asymptotic Wald type CI's are estimated for non-\code{NULL} signficance level \code{alpha}
#' for the shape parameter, assuming exactly GPD tails. When plotting on the tail index scale,
#' then a simple reciprocal transform of the CI is applied which may well be sub-optimal.
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' 
#' Pickands III, J.. (1975). Statistical inference using extreme order statistics. Annal of Statistics 3(1), 119-131.
#' 
#' Dekkers A. and de Haan, S. (1989). On the estimation of the extreme-value index and large quantile estimation.
#' Annals of Statistics 17(4), 1795-1832.
#' 
#' Resnick, S. (2007). Heavy-Tail Phenomena - Probabilistic and Statistical Modeling. Springer.
#' 
#' @author Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @section Acknowledgments: Thanks to Younes Mouatasim, Risk Dynamics, Brussels for reporting various bugs in these functions.
#' 
#' @seealso \code{\link[smoothtail:pickands]{pickands}}
#' 
#' @examples
#' \dontrun{
#' par(mfrow = c(2, 1))
#' 
#' # Reproduce graphs from Figure 4.7 of Resnick (2007)
#' data(danish, package="evir")
#' 
#' # Pickand's plot
#' pickandsplot(danish, orderlim=c(1, 150), ylim=c(-0.1, 2.2),
#'  try.thresh=c(), alpha=NULL, legend.loc=NULL)
#'  
#' # Using default settings
#' pickandsplot(danish)
#' }

pickandsplot <- function(data, orderlim = NULL, tlim = NULL,
  y.alpha = FALSE, alpha = 0.05, ylim = NULL, legend.loc = "topright",
  try.thresh = quantile(data, 0.9, na.rm = TRUE), main = "Pickand's Plot",
  xlab = "order", ylab = ifelse(y.alpha, " tail index - alpha", "shape  - xi"), ...) {

  # make sure defaults which result from function evaluations are obtained
  invisible(orderlim)
  invisible(try.thresh)
  
  # Check properties of inputs
  check.quant(data, allowna = TRUE)

  if (any(!is.finite(data))) warning("non-finite data values have been removed")
  
  # remove missing values and sort into descending order if needed
  data = data[which(is.finite(data))]
  if (is.unsorted(data)) {
    data = sort(data, decreasing = TRUE)
  } else {
    if (data[1] < data[length(data)])
      data = rev(data)
  }
  check.quant(data)

  n = length(data)
  
  # Check threshold limits if provided and set order limits if not provided
  check.param(tlim, allowvec = TRUE, allownull = TRUE)
  if (!is.null(tlim)) {
    if (length(tlim) != 2)
      stop("threshold range tlim must be a numeric vector of length 2")
    
    if (tlim[2] <= tlim[1])
      stop("a range of thresholds must be specified by tlim")

    if (is.null(orderlim)) {
      orderlim = c(sum(data >= tlim[2]), max(sum(data >= tlim[1]), 1))
    }    
  }

  # Check threshold limits if provided and set order limits if not provided
  if (!is.null(orderlim)) {
    if (length(orderlim) != 2 | mode(orderlim) != "numeric")
      stop("order statistic range orderlim must be an integer vector of length 2")

    check.n(orderlim[1])
    check.n(orderlim[2])
    
    if (orderlim[2] <= orderlim[1])
      stop("a range of order statistics must be specified by orderlim")
  
    if (orderlim[2] > floor(n/4))
      stop("maximum order statistic in orderlim must be less than floor(n/4)")
  } else {
    orderlim = c(3, floor(n/4))    
  }
  
  check.logic(y.alpha)
  check.prob(alpha, allownull = TRUE)
  if (!is.null(alpha)){
    if (alpha <= 0 | alpha >= 1)
      stop("significance level alpha must be between (0, 1)")
  }
  
  check.param(ylim, allowvec = TRUE, allownull = TRUE)
  if (!is.null(ylim)) {
    if (length(ylim) != 2) 
      stop("ylim must be a numeric vector of length 2")

    if (ylim[2] <= ylim[1])
      stop("a range of y axis limits must be specified by ylim")
  }
  
  check.text(legend.loc, allownull = TRUE)
  if (!is.null(legend.loc)) {
    if (!(legend.loc %in% c("bottomright", "bottom", "bottomleft", "left",
      "topleft", "top", "topright", "right", "center")))
      stop("legend location not correct, see help(legend)")
  }
    
  # Check given order statistics
  if (max(orderlim) <= 10)
    stop("must have more than 10 order statistics")
  
  norder = (diff(orderlim) + 1)
  if (norder < 2)
    stop("must be more than 2 order statistics considered")
  
  check.posparam(try.thresh, allowvec = TRUE, allownull = TRUE)
  if (!is.null(try.thresh)) {
    if (any((try.thresh > data[orderlim[1]]) | (try.thresh < data[orderlim[2]]))) {
      warning("potential thresholds must be within range specifed by orderlim, those outside have been set to limits")
      if (any(try.thresh > data[orderlim[1]])) {
        try.thresh[try.thresh > data[orderlim[1]]] = data[orderlim[1]]
      }
      if (any(try.thresh < data[orderlim[2]])) {
        try.thresh[try.thresh < data[orderlim[2]]] = data[orderlim[2]]
      }
      try.thresh = as.vector(try.thresh)
    }
  }

  # max order statistic
  maxks = floor(n/4)

  # order statistics
  ks = 1:maxks
  
  # Pickands estimator of xi
  Pick = log((data[ks] - data[2*ks])/(data[2*ks] - data[4*ks]))/log(2)

  # Reciprocal of Pickands estimator is tail index
  alphahat = 1/Pick
  
  # standard error of P
  Pickse = Pick*sqrt((2^(2*Pick + 1) + 1))/2/(2^Pick - 1)/log(2)/sqrt(ks)

  pickresults = data.frame(data[ks], ks, Pick, se.H = Pickse)
  
  if (!is.null(alpha)) {
    # 100(1-alpha)% CI for H
    Pickci = cbind(Pick - qnorm(1 - alpha/2) * Pickse,
                Pick + qnorm(1 - alpha/2) * Pickse)
    pickresults = cbind(pickresults, cil.Pick = Pickci[, 1], ciu.Pick = Pickci[, 2])
  }
    
  orderstats = orderlim[1]:orderlim[2]
  
  # Resolve results to be plotted
  x = ks[orderstats]
  y = Pick[orderstats]
  if (!is.null(alpha)) yci = Pickci[orderstats, ]
  norder = length(y)
  
  # xi or alpha on y-axis
  if (y.alpha) {
    y = 1/y
    if (!is.null(alpha)) yci = 1/yci
  }
  
  # Work out y-axis range (10% beyond each furthest extent of CI's)
  if (is.null(ylim)) {
    if (!is.null(alpha)) {
      ylim = range(yci, na.rm = TRUE)
    } else {
      ylim = range(y, na.rm = TRUE)      
    }
    ylim = ylim + c(-1, 1) * diff(ylim)/10
  }
  
  # Pickands plot
  par(mar = c(5, 4, 7, 2) + 0.1)
  plot(x, y, type = "l", xlab = xlab, ylab = ylab, main = main, axes = FALSE, ylim = ylim, ...)
  if (!is.null(alpha)) {
    lines(x, yci[, 1], type = "l", lty = 3)
    lines(x, yci[, 2], type = "l", lty = 3)
  }
  box()
  axis(2)

  kticks = pretty(x, 5)
  kticks = ifelse(kticks == 0, 1, kticks)
  xticks = format(data[kticks], digits = 3)

  axis(1, at = kticks)
  axis(3, at = kticks, line = 0, labels = xticks)
  mtext("Threshold", side = 3, line = 2)
  
  if (!is.null(try.thresh)) {
    ntry = length(try.thresh)
    Pickparams = rep(NA, ntry)
    linecols = rep(c("blue", "green", "red"), length.out = ntry)
    for (i in 1:ntry) {
      try.order = sum(data > try.thresh[i])      
      try.x = c(orderlim[1], try.order, orderlim[2])
      Pickparams[i] = Pick[try.order]
      try.y = rep(Pickparams[i], 3)
      
      # xi or alpha on y-axis
      if (y.alpha) try.y = 1/try.y      
      
      # Suppose to be constant above suitable threshold, different line type before and after
      lines(try.x[1:2], try.y[1:2], lwd = 2, lty = 1, col = linecols[i])
      lines(try.x[2:3], try.y[2:3], lwd = 2, lty = 2, col = linecols[i])
      abline(v = try.order, lty = 3, col = linecols[i])
    }
    if (!is.null(legend.loc)) {
      if (!is.null(alpha)) {
        legend(legend.loc, c("Pickand's Estimator", paste(100*(1 - alpha), "% CI"),
          paste("u =", formatC(try.thresh[1:min(c(3, ntry))], digits = 2, format = "g"),
            "alpha =", formatC(1/Pickparams[1:min(c(3, ntry))], digits = 2, format = "g"),
            "xi =", formatC(Pickparams[1:min(c(3, ntry))], digits = 2, format = "g"))),
          lty = c(1, 2, rep(1, min(c(3, ntry)))),
          lwd = c(2, 1, rep(1, min(c(3, ntry)))),
          col = c("black", "black", linecols), bg = "white")
      } else {
        legend(legend.loc, c("Pickand's Estimator",
          paste("u =", formatC(try.thresh[1:min(c(3, ntry))], digits = 2, format = "g"),
            "alpha =", formatC(1/Pickparams[1:min(c(3, ntry))], digits = 2, format = "g"),
            "xi =", formatC(Pickparams[1:min(c(3, ntry))], digits = 2, format = "g"))),
               lty = c(1, rep(1, min(c(3, ntry)))), lwd = c(2, rep(1, min(c(3, ntry)))),
               col = c("black", linecols), bg = "white")
      }
    }
  } else {
    if (!is.null(legend.loc)) {
      if (!is.null(alpha)) {
        legend(legend.loc, c("Pickand's Estimator", paste(100*(1 - alpha), "% CI")),
          lty = c(1, 2), lwd = c(2, 1), bg = "white")
      } else {
        legend(legend.loc, "Pickand's Estimator", lty = 1, lwd = 2, bg = "white")        
      }
    }
  }
  
  invisible(pickresults)
}
