#' @export
#' 
#' @title Hill Plot
#'
#' @description Plots the Hill plot and some its variants.
#'
#' @inheritParams mrlplot
#' @param orderlim   vector of (lower, upper) limits of order statistics
#'                   to plot estimator, or \code{NULL} to use default values
#' @param tlim       vector of (lower, upper) limits of range of threshold
#'                   to plot estimator, or \code{NULL} to use default values
#' @param hill.type  "Hill" or "SmooHill"
#' @param r          smoothing factor for "SmooHill" (integer > 1)
#' @param x.theta    logical, should order (\code{FALSE}) or theta (\code{TRUE}) be given on x-axis
#' @param y.alpha    logical, should shape xi (\code{FALSE}) or tail index alpha (\code{TRUE}) be given on y-axis
#' 
#' @details Produces the Hill, AltHill, SmooHill and AltSmooHill plots,
#'   including confidence intervals.
#'   
#'   For an ordered iid sequence \eqn{X_{(1)}\ge X_{(2)}\ge\cdots\ge X_{(n)} > 0} 
#'   the Hill (1975) estimator using \eqn{k} order statistics is given by 
#'   \deqn{H_{k,n}=\frac{1}{k}\sum_{i=1}^{k} \log(\frac{X_{(i)}}{X_{(k+1)}})}
#'   which is the pseudo-likelihood estimator of reciprocal of the tail index \eqn{\xi=/\alpha>0}
#'   for regularly varying tails (e.g. Pareto distribution).  The Hill estimator
#'   is defined on orders \eqn{k>2}, as when\eqn{k=1} the \deqn{H_{1,n}=0}. The
#'   function will calculate the Hill estimator for \eqn{k\ge 1}.
#'   The simple Hill plot is shown for \code{hill.type="Hill"}.
#'   
#'   Once a sufficiently low order statistic is reached the Hill estimator will
#'   be constant, upto sample uncertainty, for regularly varying tails. The Hill
#'   plot is a plot of \deqn{H_{k,n}} against the \eqn{k}. Symmetric asymptotic
#'   normal confidence intervals assuming Pareto tails are provided.
#'   
#'   These so called Hill's horror plots can be difficult to interpret. A smooth
#'   form of the Hill estimator was suggested by Resnick and Starica (1997): 
#'   \deqn{smooH_{k,n}=\frac{1}{(r-1)k}\sum_{j=k+1}^{rk} H_{j,n}} giving the
#'   smooHill plot which is shown for \code{hill.type="SmooHill"}. The smoothing
#'   factor is \code{r=2} by default.
#'   
#'   It has also been suggested to plot the order on a log scale, by plotting
#'   the points \eqn{(\theta, H_{\lceil n^\theta\rceil, n})} for 
#'   \eqn{0\le \theta \le 1}. This gives the so called AltHill and AltSmooHill
#'   plots. The alternative x-axis scale is chosen by \code{x.theta=TRUE}.
#'   
#'   The Hill estimator is for the GPD shape \eqn{\xi>0}, or the reciprocal of the
#'   tail index \eqn{\alpha=1/\xi>0}. The shape is plotted by default using
#'   \code{y.alpha=FALSE} and the tail index is plotted when \code{y.alpha=TRUE}.
#'   
#'   A pre-chosen threshold (or more than one) can be given in
#'   \code{try.thresh}. The estimated parameter (\eqn{\xi} or \eqn{\alpha}) at
#'   each threshold are plot by a horizontal solid line for all higher thresholds. 
#'   The threshold should be set as low as possible, so a dashed line is shown
#'   below the pre-chosen threshold. If the Hill estimator is similar to the
#'   dashed line then a lower threshold may be chosen.
#'   
#'   If no order statistic (or threshold) limits are provided \code{orderlim =
#'   tlim = NULL} then the lowest order statistic is set to \eqn{X_{(3)}} and
#'   highest possible value \eqn{X_{(n-1)}}. However, the Hill estimator is always
#'   output for all \eqn{k=1, \ldots, n-1} and \eqn{k=1, \ldots, floor(n/k)} for
#'   smooHill estimator.
#'   
#'   The missing (\code{NA} and \code{NaN}) and non-finite values are ignored.
#'   Non-positive data are ignored.
#'   
#'   The lower x-axis is the order \eqn{k} or \eqn{\theta}, chosen by the option
#'   \code{x.theta=FALSE} and \code{x.theta=TRUE} respectively. The upper axis
#'   is for the corresponding threshold.
#' 
#' @return \code{\link[evmix:hillplot]{hillplot}} gives the Hill plot. It also 
#'   returns a dataframe containing columns of the order statistics, order, Hill
#'   estimator, it's standard devation and \eqn{100(1 - \alpha)\%} confidence
#'   interval (when requested). When the SmooHill plot is selected, then the corresponding
#'   SmooHill estimates are appended.
#' 
#' @note 
#' Warning: Hill plots are not location invariant.
#' 
#' Asymptotic Wald type CI's are estimated for non-\code{NULL} signficance level \code{alpha}
#' for the shape parameter, assuming exactly Pareto tails. When plotting on the tail index scale,
#' then a simple  reciprocal transform of the CI is applied which may be sub-optimal.
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' 
#' Hill, B.M. (1975). A simple general approach to inference about the tail of a distribution. Annals of Statistics 13, 331-341.
#' 
#' Resnick, S. and Starica, C. (1997). Smoothing the Hill estimator. Advances in Applied Probability 29, 271-293.
#' 
#' Resnick, S. (1997). Discussion of the Danish Data of Large Fire Insurance Losses. Astin Bulletin 27, 139-151.
#' 
#' @author Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @section Acknowledgments: Thanks to Younes Mouatasim, Risk Dynamics, Brussels for reporting various bugs in these functions.
#' 
#' @seealso \code{\link[evir:hill]{hill}}
#' 
#' @examples
#' \dontrun{
#' # Reproduce graphs from Figure 2.4 of Resnick (1997)
#' data(danish, package="evir")
#' par(mfrow = c(2, 2))
#' 
#' # Hill plot
#' hillplot(danish, y.alpha=TRUE, ylim=c(1.1, 2))
#' 
#' # AltHill plot
#' hillplot(danish, y.alpha=TRUE, x.theta=TRUE, ylim=c(1.1, 2))
#' 
#' # AltSmooHill plot
#' hillplot(danish, hill.type="SmooHill", r=3, y.alpha=TRUE, x.theta=TRUE, ylim=c(1.35, 1.85))
#' 
#' # AltHill and AltSmooHill plot (no CI's or legend)
#' hillout = hillplot(danish, hill.type="SmooHill", r=3, y.alpha=TRUE, 
#'  x.theta=TRUE, try.thresh = c(), alpha=NULL, ylim=c(1.1, 2), legend.loc=NULL, lty=2)
#' n = length(danish)
#' with(hillout[3:n,], lines(log(ks)/log(n), 1/H, type="s"))
#' }

hillplot <- function(data, orderlim = NULL, tlim = NULL, hill.type = "Hill", r = 2,
  x.theta = FALSE, y.alpha = FALSE, alpha = 0.05, ylim = NULL, legend.loc = "topright",
  try.thresh = quantile(data[data > 0], 0.9, na.rm = TRUE),
  main = paste(ifelse(x.theta, "Alt", ""), hill.type, " Plot", sep=""),
  xlab = ifelse(x.theta, "theta", "order"),
  ylab = paste(ifelse(x.theta, "Alt", ""), hill.type, ifelse(y.alpha, " alpha", " xi"), ">0", sep=""),
  ...) {

  # make sure defaults which result from function evaluations are obtained
  invisible(orderlim)
  invisible(try.thresh)
  
  # Check properties of inputs
  check.quant(data, allowna = TRUE, allowinf = TRUE)

  if (any(!is.finite(data))) warning("non-finite have been removed")
  
  # remove missing values
  data = data[which(is.finite(data))]

  if (any(data <= 0)) {
    warning("non-positive values have been removed")
    data = data[data > 0]
  }
  check.quant(data)

  # sort into descending order if needed
  if (is.unsorted(data)) {
    data = sort(data, decreasing = TRUE)
  } else {
    if (data[1] < data[length(data)])
      data = rev(data)
  }
  n = length(data)
  
  check.text(hill.type)
  if (!(hill.type %in% c("Hill", "SmooHill")))
    stop("Hill plot type must be Hill or SmooHill, see help(hillplot)")
  
  # check r if smoothing Hill plot
  if (hill.type == "SmooHill") {
    check.n(r)
    
    if (r >= (n/10)) stop("smoothing factor r must be an integer > 1, typically 2 or 3")
  }

  # Check threshold limits if provided and set order limits if not provided
  check.posparam(tlim, allowvec = TRUE, allownull = TRUE)
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
    if (length(orderlim) != 2)
      stop("order statistic range orderlim must be an integer vector of length 2")

    check.n(orderlim[1])
    check.n(orderlim[2])
    
    if (orderlim[2] <= orderlim[1])
      stop("a range of order statistics must be specified by orderlim")
  
    if ((hill.type == "Hill") & (orderlim[2] > (n - 1)))
      stop("maximum order statistic in orderlim must be less then n-1")
    
    if ((hill.type == "SmooHill") & (orderlim[2] > floor(n/r)))
      stop("maximum order statistic in orderlim must be less then floor(n/r)")
  } else {
    orderlim = c(3, ifelse(hill.type == "SmooHill", floor(n/r), n - 1))
  }

  check.logic(x.theta)
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
    stop("must plot more than 10 order statistics")
  
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
  
  logdata = log(data)
  
  # order statistics
  ks = 1:(n-1)
  
  # logarithm scale for order statistic
  theta = log(ks)/log(n)

  # Hill estimator of xi
  H = cumsum(logdata[-n])/ks - logdata[-1]

  # Reciprocal of Hill estimator is tail index
  alphahat = 1/H
  
  # standard error of H
  Hse = H/sqrt(ks)

  hillresults = data.frame(data[-n], ks, H, se.H = Hse)
  
  if (!is.null(alpha)) {
    # 100(1-alpha)% Wald CI for H
    Hci = cbind(H - qnorm(1 - alpha/2) * Hse,
                H + qnorm(1 - alpha/2) * Hse)
    hillresults = cbind(hillresults, cil.H = Hci[, 1], ciu.H = Hci[, 2])
  }

  # smooHill estimator, if needed
  if (hill.type == "SmooHill") {
    # max order statistic for given smoothing factor
    smook = floor(n/r)
    
    # order statistics and smooHill estimator of xi
    smooks = 1:smook
    Hc = cumsum(H)
    smooH = (Hc[r*smooks] - Hc[smooks])/(r - 1)/smooks
  
    # standard error of Hsmoo
    smooHse = smooH * sqrt(2*(1 - (log(r)/(r - 1)))/(r - 1))/sqrt(smooks)

    hillresults = cbind(hillresults, smooks = NA, smooH = NA, se.smooH = NA)
    hillresults$smooks[smooks] = smooks
    hillresults$smooH[smooks] = smooH
    hillresults$se.smooH[smooks] = smooHse
    
    if (!is.null(alpha)) {
      # 100(1-alpha)% CI for Hsmoo
      smooHci = cbind(smooH - qnorm(1 - alpha/2) * smooHse,
                  smooH + qnorm(1 - alpha/2) * smooHse)
      hillresults = cbind(hillresults, cil.smooH = NA, ciu.smooH = NA)
      hillresults$cil.smooH[smooks] = smooHci[, 1]
      hillresults$ciu.smooH[smooks] = smooHci[, 2]
    }        
  }

  # Resolve results to be plotted
  if (hill.type == "Hill") {
    x = ks
    y = H
    if (!is.null(alpha)) yci = Hci
  } else {
    x = smooks
    y = smooH
    if (!is.null(alpha)) yci = smooHci    
  }
  orderstats = orderlim[1]:orderlim[2]
  
  x = x[orderstats]
  y = y[orderstats]
  if (!is.null(alpha)) yci = yci[orderstats, ]
  norder = length(y)
  
  # xi or alpha on y-axis
  if (y.alpha) {
    y = 1/y
    if (!is.null(alpha)) yci = 1/yci # sub-optimal?
  }
  
  # order statistic or theta on the usual x-axis
  if (x.theta) x = log(x)/log(n)
  
  # Work out y-axis range (10% beyond each furthest extent of CI's)
  if (is.null(ylim)) {
    if (!is.null(alpha)) {
      ylim = range(yci, na.rm = TRUE)
    } else {
      ylim = range(y, na.rm = TRUE)
    }
    ylim = ylim + c(-1, 1) * diff(ylim)/10
  }
  
  # Hill plot (or it's variants)
  par(mar = c(5, 4, 7, 2) + 0.1)
  plot(x, y, type = ifelse(x.theta, "s", "l"),
       xlab = xlab, ylab = ylab, main = main, axes = FALSE, ylim = ylim, ...)
  if (!is.null(alpha)) {
    lines(x, yci[, 1], type = ifelse(x.theta, "s", "l"), lty = 3)
    lines(x, yci[, 2], type = ifelse(x.theta, "s", "l"), lty = 3)
  }
  box()
  axis(2)

  kticks = pretty(x, 5)
  if (x.theta) {
    xticks = prettyNum(format(data[floor(n^kticks)],digits = 3), drop0trailing = TRUE)
  } else {
    kticks = ifelse(kticks == 0, 1, kticks)
    xticks = format(data[kticks], digits = 3)
  }
  axis(1, at = kticks)
  axis(3, at = kticks, line = 0, labels = xticks)
  mtext("Threshold", side = 3, line = 2)
  
  legend.text = paste(ifelse(x.theta, "Alt", ""), hill.type, sep="")
  
  if (!is.null(try.thresh)) {
    ntry = length(try.thresh)
    Hparams = rep(NA, ntry)
    linecols = rep(c("blue", "green", "red"), length.out = ntry)
    for (i in 1:ntry) {
      try.order = sum(data > try.thresh[i])
      try.x = c(orderlim[1], try.order, orderlim[2])
      Hparams[i] = ifelse(hill.type == "SmooHill", smooH[try.order], H[try.order])
      try.y = rep(Hparams[i], 3)
      
      # xi or alpha on y-axis
      if (y.alpha) try.y = 1/try.y      
      
      # order statistic or theta on the usual x-axis
      if (x.theta) {
        try.x = log(try.x)/log(n)
        try.order = log(try.order)/log(n)
      }

      # Suppose to be constant above suitable threshold, different line type before and after
      lines(try.x[1:2], try.y[1:2], lwd = 2, lty = 1, col = linecols[i])
      lines(try.x[2:3], try.y[2:3], lwd = 2, lty = 2, col = linecols[i])
      abline(v = try.order, lty = 3, col = linecols[i])
    }
    if (!is.null(legend.loc)) {
      if (!is.null(alpha)) {
        legend(legend.loc, c(legend.text, paste(100*(1 - alpha), "% CI"),
          paste("u =", formatC(try.thresh[1:min(c(3, ntry))], digits = 2, format = "g"),
            "alpha =", formatC(1/Hparams[1:min(c(3, ntry))], digits = 2, format = "g"),
            "xi =", formatC(Hparams[1:min(c(3, ntry))], digits = 2, format = "g"))),
          lty = c(1, 2, rep(1, min(c(3, ntry)))),
          lwd = c(2, 1, rep(1, min(c(3, ntry)))),
          col = c("black", "black", linecols), bg = "white")
      } else {
        legend(legend.loc, c(legend.text,
          paste("u =", formatC(try.thresh[1:min(c(3, ntry))], digits = 2, format = "g"),
            "alpha =", formatC(1/Hparams[1:min(c(3, ntry))], digits = 2, format = "g"),
            "xi =", formatC(Hparams[1:min(c(3, ntry))], digits = 2, format = "g"))),
               lty = c(1, rep(1, min(c(3, ntry)))), lwd = c(2, rep(1, min(c(3, ntry)))),
               col = c("black", linecols), bg = "white")
      }
    }
  } else {
    if (!is.null(legend.loc)) {
      if (!is.null(alpha)) {
        legend(legend.loc, c(legend.text, paste(100*(1 - alpha), "% CI")),
          lty = c(1, 2), lwd = c(2, 1), bg = "white")
      } else {
        legend(legend.loc, legend.text, lty = 1, lwd = 2, bg = "white")        
      }
    }
  }
  
  invisible(hillresults)
}
