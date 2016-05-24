#' @export
#' 
#' @title Diagnostic Plots for Extreme Value Mixture Models
#'
#' @description The classic four diagnostic plots for evaluating extreme value mixture models: 
#'   1) return level plot, 2) Q-Q plot, 3) P-P plot and 4) density plot. Each plot is available
#'   individually or as the usual 2x2 collection.
#'
#' @param modelfit   fitted extreme value mixture model object
#' @param upperfocus logical, should plot focus on upper tail?
#' @param alpha      significance level over range (0, 1), or \code{NULL} for no CI
#' @param N          number of Monte Carlo simulation for CI (N>=10)
#' @param legend     logical, should legend be included
#' @param ...        further arguments to be passed to the plotting functions
#' 
#' @details Model diagnostics are available for all the fitted extreme mixture models in the 
#' \code{\link[evmix:evmix-package]{evmix}} package. These \code{modelfit} is output by all the fitting 
#' functions, e.g. \code{\link[evmix:fgpd]{fgpd}} and \code{\link[evmix:fnormgpd]{fnormgpd}}.
#' 
#' Consistent with \code{\link[evd:plot.uvevd]{plot}} function in the 
#' \code{\link[evd:gpd]{evd}} library the \code{\link[stats:ppoints]{ppoints}} to 
#' estimate the empirical cumulative probabilities. The default behaviour of this
#' function is to use \deqn{(i-0.5)/n} as the estimate for the \eqn{i}th order statistic of
#' the given sample of size \eqn{n}.
#' 
#' The return level plot has the quantile (\eqn{q} where \eqn{P(X \ge q)=p} on
#' the \eqn{y}-axis, for a particular survival probability \eqn{p}. The return period
#' \eqn{t=1/p} is shown on the \eqn{x}-axis. The return level is given by:
#'    \deqn{q = u + \sigma_u [(\phi_u t)^\xi - 1]/\xi}
#' for \eqn{\xi\ne 0}. But in the case of \eqn{\xi = 0} this simplifies to 
#'  \deqn{q = u + \sigma_u log(\phi_u t)}
#' which is linear when plotted against the return period on a logarithmic scale. The special
#' case of exponential/Type I (\eqn{\xi=0}) upper tail behaviour will be linear on
#' this scale. This is the same tranformation as in the GPD/POT diagnostic plot function
#' \code{\link[evd:plot.uvevd]{plot.uvevd}} in the \code{\link[evd:plot.uvevd]{evd}} package,
#' from which these functions were derived.
#' 
#' The crosses are the empirical quantiles/return levels (i.e. the ordered sample data)
#' against their corresponding transformed empirical return period (from 
#' \code{\link[stats:ppoints]{ppoints}}). The solid line is the theoretical return level
#' (quantile) function using the estimated parameters. The estimated threshold 
#' \code{u} and tail fraction \code{phiu} are shown. For the two tailed models both
#' thresholds \code{ul} and \code{ur} and corresponding tail fractions 
#' \code{phiul} and \code{phiur} are shown. The approximate pointwise confidence intervals
#' for the quantiles are obtained by Monte Carlo simulation using the estimated parameters.
#' Notice that these intervals ignore the parameter estimation uncertainty.
#' 
#' The Q-Q and P-P plots have the empirical values on the \eqn{y}-axis and theoretical values
#' from the fitted model on the \eqn{x}-axis.
#' 
#' The density plot provides a histogram of the sample data overlaid with the fitted density
#' and a standard kernel density estimate using the \code{\link[stats:density]{density}}
#' function. The default settings for the \code{\link[stats:density]{density}} function are used.
#' Note that for distributions with bounded support (e.g. GPD) with high density near the
#' boundary standard kernel density estimators exhibit a negative bias due to leakage past
#' the boundary. So in this case they should not be taken too seriously.
#' 
#' For the kernel density estimates (i.e. \code{kden} and \code{bckden}) there is no threshold, 
#' so no upper tail focus is carried out.
#' 
#' See \code{\link[evd:plot.uvevd]{plot.uvevd}} for more detailed explanations of these
#' types of plots.
#' 
#' @return \code{\link[evmix:evmix.diag]{rlplot}} gives the return level plot, 
#' \code{\link[evmix:evmix.diag]{qplot}} gives the Q-Q plot,
#' \code{\link[evmix:evmix.diag]{pplot}} gives the P-P plot,
#' \code{\link[evmix:evmix.diag]{densplot}} gives density plot and
#' \code{\link[evmix:evmix.diag]{evmix.diag}} gives the collection of all 4.
#' 
#' @note For all mixture models the missing values are removed by the fitting functions 
#' (e.g. \code{\link[evmix:fnormgpd]{fnormgpd}} and \code{\link[evmix:fgng]{fgng}}).
#' However, these are retained in the GPD fitting \code{\link[evmix:fgpd]{fgpd}}, as they 
#' are interpreted as values below the threshold.
#' 
#' By default all the plots focus in on the upper tail, but they can be used to display 
#' the fit over the entire range of support.
#' 
#' You cannot pass \code{xlim} or \code{ylim} to the plotting functions via \code{...}
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' 
#' \url{http://en.wikipedia.org/wiki/Q-Q_plot}
#' 
#' \url{http://en.wikipedia.org/wiki/P-P_plot}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Coles S.G. (2004). An Introduction to the Statistical Modelling of Extreme Values.
#' Springer-Verlag: London.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @section Acknowledgments: Based on the GPD/POT diagnostic function \code{\link[evd:plot.uvevd]{plot.uvevd}} in the \code{\link[evd:plot.uvevd]{evd}} package for which Stuart Coles' and Alec Stephenson's 
#' contributions are gratefully acknowledged.
#' They are designed to have similar syntax and functionality to simplify the transition for users of these packages.
#' 
#' @seealso \code{\link[stats:ppoints]{ppoints}}, \code{\link[evd:plot.uvevd]{plot.uvevd}} and
#'   \code{\link[ismev:gpd.diag]{gpd.diag}}.
#' @aliases  evmix.diag rlplot qplot pplot densplot
#' @family   evmix.diag
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' 
#' x = sort(rnorm(1000))
#' fit = fnormgpd(x)
#' evmix.diag(fit)
#' 
#' # repeat without focussing on upper tail
#' par(mfrow=c(2,2))
#' rlplot(fit, upperfocus = FALSE)
#' qplot(fit, upperfocus = FALSE)
#' pplot(fit, upperfocus = FALSE)
#' densplot(fit, upperfocus = FALSE)
#' }
#' 

evmix.diag <- function(modelfit, upperfocus = TRUE, alpha = 0.05, N = 1000,
  legend = FALSE, ...) {
  
  # Check properties of inputs
  if (missing(modelfit))
    stop("modelfit from extreme value mixture model must be given")
   
  if (!exists("call", modelfit))
    stop("modelfit from extreme value mixture model must be given")

  modelname = substring(strsplit(deparse(modelfit$call),"\\(")[[c(1, 1)]], 2)

  allmodelnames = c("gpd", "betagpd", "betagpdcon", 
    "bckden", "bckdengpd", "bckdengpdcon", "dwm", "gammagpd", "gammagpdcon",
    "gkg", "gkgcon", "gng", "gngcon", "hpd", "hpdcon", 
    "itmnormgpd", "itmweibullgpd", "itmgng", "mgammagpd", 
    "kden", "kdengpd", "kdengpdcon", "lognormgpd", "lognormgpdcon",
    "normgpd", "normgpdcon", "weibullgpd", "weibullgpdcon")

  if (!(modelname %in% allmodelnames))
    stop("invalid extreme value mixture model given")

  check.logic(upperfocus)
  check.logic(legend)
  check.prob(alpha, allownull = TRUE)
  if (!is.null(alpha)) 
    if ((alpha == 0) | (alpha == 1)) stop("alpha must be in (0, 1)")
  check.n(N)

  if (is.unsorted(modelfit$x, na.rm = TRUE)) {
    modelfit$x=sort(modelfit$x)
  } else {
    if (modelfit$x[1] > modelfit$x[length(modelfit$x)])
      modelfit$x = rev(modelfit$x)
  }
  
  par(mfrow = c(2, 2))
  rlplot(modelfit, upperfocus, alpha, N, legend, ...)
  qplot(modelfit, upperfocus, alpha, N, legend, ...)
  pplot(modelfit, upperfocus, alpha, N, legend, ...)
  densplot(modelfit, upperfocus, legend, ...)
}

#' @export
#' @aliases evmix.diag rlplot qplot pplot densplot
#' @rdname evmix.diag

rlplot <- function(modelfit, upperfocus = TRUE, alpha = 0.05, N = 1000,
  legend = TRUE, ...) {
  
  # Check properties of inputs
  if (missing(modelfit))
    stop("modelfit from extreme value mixture model must be given")
   
  if (!exists("call", modelfit))
    stop("modelfit from extreme value mixture model must be given")

  modelname = substring(strsplit(deparse(modelfit$call),"\\(")[[c(1, 1)]], 2)

  allmodelnames = c("gpd", "betagpd", "betagpdcon", 
    "bckden", "bckdengpd", "bckdengpdcon", "dwm", "gammagpd", "gammagpdcon",
    "gkg", "gkgcon", "gng", "gngcon", "hpd", "hpdcon", 
    "itmnormgpd", "itmweibullgpd", "itmgng", "mgammagpd", 
    "kden", "kdengpd", "kdengpdcon", "lognormgpd", "lognormgpdcon",
    "normgpd", "normgpdcon", "psden", "psdengpd", "weibullgpd", "weibullgpdcon")

  if (!(modelname %in% allmodelnames))
    stop("invalid extreme value mixture model given")

  check.logic(upperfocus)
  check.logic(legend)
  check.prob(alpha, allownull = TRUE)
  if (!is.null(alpha)) 
    if ((alpha == 0) | (alpha == 1)) stop("alpha must be in (0, 1)")
  check.n(N)

  if (is.unsorted(modelfit$x, na.rm = TRUE)) {
    modelfit$x=sort(modelfit$x)
  } else {
    if (modelfit$x[1] > modelfit$x[length(modelfit$x)])
      modelfit$x = rev(modelfit$x)
  }
  
  nothresh = (modelname == "kden") | (modelname =="bckden") | (modelname =="dwm") | (modelname == "psden")
  twotail = (modelname == "gng") | (modelname == "gngcon") | (modelname == "gkg") | (modelname == "gkgcon")
  nophiu = (modelname == "itmnormgpd") | (modelname == "itmweibullgpd") | (modelname == "itmgng")
  
  if (nothresh) {
    upperfocus = FALSE
  }
  if (modelname == "gpd") {
    upperfocus = TRUE  
  }

  # makes following code easier to read
  if (twotail) {
    phiu = modelfit$phiur
    u = modelfit$ur
  } else if (nothresh) {
    phiu = -Inf
    u = Inf 
  } else if (nophiu) {
    phiu = -Inf
    u = ifelse(modelname == "itmgng", modelfit$ur, modelfit$u)
  } else {
    phiu = modelfit$phiu
    u = modelfit$u    
  }
  
  n = length(modelfit$x)
  emp.prob = ppoints(n)
  trans.emp.prob = 1/(1 - emp.prob)

  min.emp.power = -log10(1 - 1/n/100)
  max.emp.power = ceiling(log10(max(trans.emp.prob))) + 1

  the.prob = 1 - 10^(-seq(min.emp.power, max.emp.power, 0.05))
  # Need both tail to be in fine detail
  the.prob = sort(c(the.prob, 1 - the.prob)) 
  trans.the.prob = 1/(1-the.prob)

  the.quant = switch(modelname,
    gpd = qgpd(the.prob, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    bckden = qbckden(the.prob, modelfit$kerncentres, modelfit$lambda, bw = NULL,
      modelfit$kernel, modelfit$bcmethod,
      modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
    bckdengpd = qbckdengpd(the.prob, modelfit$x, modelfit$lambda, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel,
      modelfit$bcmethod, modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
    bckdengpdcon = qbckdengpdcon(the.prob, modelfit$x, modelfit$lambda, modelfit$u,
      modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel,
      modelfit$bcmethod, modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
    betagpd = qbetagpd(the.prob, modelfit$bshape1, modelfit$bshape2, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    betagpdcon = qbetagpdcon(the.prob, modelfit$bshape1, modelfit$bshape2, modelfit$u,
      modelfit$xi, modelfit$phiu),
    dwm = qdwm(the.prob, modelfit$wshape, modelfit$wscale, modelfit$cmu, modelfit$ctau,
      modelfit$sigmau, modelfit$xi),
    gammagpd = qgammagpd(the.prob, modelfit$gshape, modelfit$gscale, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    gammagpdcon = qgammagpdcon(the.prob, modelfit$gshape, modelfit$gscale, modelfit$u,
      modelfit$xi, modelfit$phiu),
    gkg = qgkg(the.prob, modelfit$x, modelfit$lambda, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul,
      modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur, bw = NULL, modelfit$kernel),
    gkgcon = qgkgcon(the.prob, modelfit$x, modelfit$lambda, 
      modelfit$ul, modelfit$xil, modelfit$phiul,
      modelfit$ur, modelfit$xir, modelfit$phiur, bw = NULL, modelfit$kernel),
    gng = qgng(the.prob, modelfit$nmean, modelfit$nsd, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
      modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
    gngcon = qgngcon(the.prob, modelfit$nmean, modelfit$nsd, 
      modelfit$ul, modelfit$xil, modelfit$phiul, modelfit$ur, modelfit$xir, modelfit$phiur),
    hpd = qhpd(the.prob, modelfit$nmean, modelfit$nsd, modelfit$xi),
    hpdcon = qhpdcon(the.prob, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$xi),
    itmnormgpd = qitmnormgpd(the.prob, modelfit$nmean, modelfit$nsd, modelfit$epsilon, 
      modelfit$u, modelfit$sigmau, modelfit$xi),
    itmweibullgpd = qitmweibullgpd(the.prob, modelfit$wshape, modelfit$wscale, modelfit$epsilon, 
      modelfit$u, modelfit$sigmau, modelfit$xi),
    itmgng = qitmgng(the.prob, modelfit$nmean, modelfit$nsd, modelfit$epsilon, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$ur, modelfit$sigmaur, modelfit$xir),
    kden = qkden(the.prob, modelfit$kerncentres, modelfit$lambda, bw = NULL, modelfit$kernel),
    kdengpd = qkdengpd(the.prob, modelfit$x, modelfit$lambda, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel),
    kdengpdcon = qkdengpdcon(the.prob, modelfit$x, modelfit$lambda, modelfit$u,
      modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel),
    mgammagpd = qmgammagpd(the.prob, modelfit$mgshape, modelfit$mgscale, modelfit$mgweight,
      modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    lognormgpd = qlognormgpd(the.prob, modelfit$lnmean, modelfit$lnsd, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    lognormgpdcon = qlognormgpd(the.prob, modelfit$lnmean, modelfit$lnsd, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu), 
    normgpd = qnormgpd(the.prob, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau,
      modelfit$xi, modelfit$phiu),
    normgpdcon = qnormgpdcon(the.prob, modelfit$nmean, modelfit$nsd, modelfit$u,
      modelfit$xi, modelfit$phiu),
    psden = qpsden(the.prob, modelfit$beta, modelfit$nbinwidth, modelfit$xrange, modelfit$nseg,
                   modelfit$degree, modelfit$design.knots),
    psdengpd = qpsdengpd(the.prob, modelfit$beta, modelfit$nbinwidth, modelfit$xrange, modelfit$nseg,
                   modelfit$degree, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu, modelfit$design.knots),
    weibullgpd = qweibullgpd(the.prob, modelfit$wshape, modelfit$wscale, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    weibullgpdcon = qweibullgpdcon(the.prob, modelfit$wshape, modelfit$wscale, modelfit$u,
      modelfit$xi, modelfit$phiu))

  # If focus on upper tail then must display at least upper 10%
  # even if tail fraction is smaller, otherwise may not look nice
  if (upperfocus) { 
    xlims = c(1/ifelse(phiu < 0.1, 0.1, phiu), 10^max.emp.power)
    ylims = c(min(quantile(modelfit$x, 0.9, na.rm = TRUE), u), 
      max(c(modelfit$x, the.quant), na.rm = TRUE))
  } else {
    xlims = c(min(trans.emp.prob), 10^max.emp.power)
    ylims = range(c(modelfit$x, the.quant), na.rm = TRUE)
  }

  if (!is.null(alpha)) {
    # obtain CI by Monte Carlo simulation (ignores estimation uncertainty)
    simulate.new.quantiles <- function(i, modelfit, modelname) {
      simdata = switch(modelname,
        gpd = rgpd(n, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        bckden = rbckden(n, modelfit$kerncentres, modelfit$lambda, bw = NULL, modelfit$kernel,
          modelfit$bcmethod, modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
        bckdengpd = rbckdengpd(n, modelfit$x, modelfit$lambda, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel,
          modelfit$bcmethod, modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
        bckdengpdcon = rbckdengpdcon(n, modelfit$x, modelfit$lambda, modelfit$u,
          modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel,
          modelfit$bcmethod, modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
        betagpd = rbetagpd(n, modelfit$bshape1, modelfit$bshape2, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        betagpdcon = rbetagpdcon(n, modelfit$bshape1, modelfit$bshape2, modelfit$u,
          modelfit$xi, modelfit$phiu),
        dwm = rdwm(n, modelfit$wshape, modelfit$wscale, modelfit$cmu, modelfit$ctau,
          modelfit$sigmau, modelfit$xi),
        gammagpd = rgammagpd(n, modelfit$gshape, modelfit$gscale, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        gammagpdcon = rgammagpdcon(n, modelfit$gshape, modelfit$gscale, modelfit$u,
          modelfit$xi, modelfit$phiu),
        gkg = rgkg(n, modelfit$x, modelfit$lambda, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul,
          modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur, bw = NULL, modelfit$kernel),
        gkgcon = rgkgcon(n, modelfit$x, modelfit$lambda, 
          modelfit$ul, modelfit$xil, modelfit$phiul,
          modelfit$ur, modelfit$xir, modelfit$phiur, bw = NULL, modelfit$kernel),
        gng = rgng(n, modelfit$nmean, modelfit$nsd, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
          modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
        gngcon = rgngcon(n, modelfit$nmean, modelfit$nsd, 
          modelfit$ul, modelfit$xil, modelfit$phiul,
          modelfit$ur, modelfit$xir, modelfit$phiur),
        hpd = rhpd(n, modelfit$nmean, modelfit$nsd, modelfit$xi),
        hpdcon = rhpdcon(n, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$xi),
        itmnormgpd = ritmnormgpd(n, modelfit$nmean, modelfit$nsd, modelfit$epsilon, 
          modelfit$u, modelfit$sigmau, modelfit$xi),
        itmweibullgpd = ritmweibullgpd(n, modelfit$wshape, modelfit$wscale, modelfit$epsilon, 
          modelfit$u, modelfit$sigmau, modelfit$xi),
        itmgng = ritmgng(n, modelfit$nmean, modelfit$nsd, modelfit$epsilon, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$ur, modelfit$sigmaur, modelfit$xir),
        kden = rkden(n, modelfit$kerncentres, modelfit$lambda, bw = NULL, modelfit$kernel),
        kdengpd = rkdengpd(n, modelfit$x, modelfit$lambda, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel),
        kdengpdcon = rkdengpdcon(n, modelfit$x, modelfit$lambda, modelfit$u,
          modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel),
        lognormgpd = rlognormgpd(n, modelfit$lnmean, modelfit$lnsd, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        lognormgpdcon = rlognormgpd(n, modelfit$lnmean, modelfit$lnsd, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        mgammagpd = rmgammagpd(n, modelfit$mgshape, modelfit$mgscale, modelfit$mgweight,
          modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        normgpd = rnormgpd(n, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau,
          modelfit$xi, modelfit$phiu),
        normgpdcon = rnormgpdcon(n, modelfit$nmean, modelfit$nsd, modelfit$u,
          modelfit$xi, modelfit$phiu),
        psden = rpsden(n, modelfit$beta, modelfit$nbinwidth, modelfit$xrange, modelfit$nseg,
                   modelfit$degree, modelfit$design.knots),
        psdengpd = rpsdengpd(n, modelfit$beta, modelfit$nbinwidth, modelfit$xrange, modelfit$nseg,
                   modelfit$degree, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu, modelfit$design.knots),
        weibullgpd = rweibullgpd(n, modelfit$wshape, modelfit$wscale, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        weibullgpdcon = rweibullgpdcon(n, modelfit$wshape, modelfit$wscale, modelfit$u,
          modelfit$xi, modelfit$phiu))
      sort(simdata, na.last = FALSE)
    }
    sim.q = sapply(1:N, FUN = simulate.new.quantiles, modelfit = modelfit, modelname = modelname)
    ci.q = apply(sim.q, FUN = quantile, MARGIN = 1, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)

    # simulation for gpd allow for uncertainty in phiu, so need to ignore CI below threshold
    if (modelname == "gpd") {
      ci.q[1, modelfit$x < u] = NA
      ci.q[2, modelfit$x < u] = NA
    }
    ylims[2] = max(c(ci.q, ylims[2]), na.rm = TRUE)
  }

  if (is.infinite(ylims[2])) ylims[2] = max(modelfit$x)

  plot(trans.emp.prob, modelfit$x, pch = 'x', cex = 0.8, log = "x",
    xlab = "Return Period", ylab = "Return Level",
    xlim = xlims, ylim = ylims, axes = FALSE, ...)
  lines(trans.the.prob, the.quant, col="red", lwd = 2)

  xpower = seq(0, max.emp.power)
  axis(1, at = 10^xpower, labels = formatC(10^xpower, format = "d"))
  axis(2)
  box()


  if ((twotail) | (modelname == "itmgng")) {
    abline(h = c(modelfit$ul, modelfit$ur), lty = 3)
  } else {
    abline(h = u, lty = 3)
  }

  if (!((nothresh) | (nophiu)))  abline(v = 1/phiu, lty = 3)
  
  if (!is.null(alpha)) {
    lines(trans.emp.prob, ci.q[1, ], lty=2)
    lines(trans.emp.prob, ci.q[2, ], lty=2)
    if (legend) {
      if (nophiu) {
        legend("bottomright", c("Data", paste("Fitted", modelname, "Model"),
          paste("Simulated Pointwise ", formatC(100 * (1 - alpha), format = "f", digits=1), "% CI", sep=""),
          ifelse((modelname == "itmgng") & !upperfocus, "Thresholds", "Threshold")),
          pch = c(4, rep(-1, 3)), lty = c(0, 1, 2, 3), bg = "white", col = c("black", "red", "black"))
      } else if (nothresh) {
        legend("bottomright", c("Data", paste("Fitted", modelname, "Model"),
          paste("Simulated Pointwise ", formatC(100 * (1 - alpha), format = "f", digits=1), "% CI", sep="")),
          pch = c(4, rep(-1, 2)), lty = c(0, 1, 2), bg = "white", col = c("black", "red"))
      } else {
        legend("bottomright", c("Data", paste("Fitted", modelname, "Model"),
          paste("Simulated Pointwise ", formatC(100 * (1 - alpha), format = "f", digits=1), "% CI", sep=""),
          ifelse(twotail & !upperfocus,
            "Tail Fractions and Thresholds", "Upper Tail Fraction and Threshold")),
          pch = c(4, rep(-1, 3)), lty = c(0, 1, 2, 3), bg = "white", col = c("black", "red", "black"))
      }
    }
  } else {
    if (legend) {
      if (nophiu) {
        legend("bottomright", c("Data", paste("Fitted", modelname, "Model"),
          ifelse((modelname == "itmgng") & !upperfocus, "Thresholds", "Threshold")),
          pch = c(4, rep(-1, 2)), lty = c(0, 1, 3), bg = "white", col = c("black", "red", "black"))
      } else if (nothresh) {
        legend("bottomright", c("Data", paste("Fitted", modelname, "Model")),
          pch = c(4, -1), lty = c(0, 1), bg = "white", col = c("black", "red"))
      } else {
        legend("bottomright", c("Data", paste("Fitted", modelname, "Model"),
          ifelse((twotail) & !upperfocus,
            "Tail Fractions and Thresholds", "Upper Tail Fraction and Threshold")),
          pch = c(4, rep(-1, 2)), lty = c(0, 1, 3), bg = "white", col = c("black", "red", "black"))
      }
    }
  }
}

#' @export
#' @aliases evmix.diag rlplot qplot pplot densplot
#' @rdname evmix.diag

qplot <- function(modelfit, upperfocus = TRUE, alpha = 0.05, N = 1000, 
  legend = TRUE, ...) {

  # Check properties of inputs
  if (missing(modelfit))
    stop("modelfit from extreme value mixture model must be given")
   
  if (!exists("call", modelfit))
    stop("modelfit from extreme value mixture model must be given")

  modelname = substring(strsplit(deparse(modelfit$call), "\\(")[[c(1, 1)]], 2)

  allmodelnames = c("gpd", "betagpd", "betagpdcon", 
    "bckden", "bckdengpd", "bckdengpdcon", "dwm", "gammagpd", "gammagpdcon",
    "gkg", "gkgcon", "gng", "gngcon", "hpd", "hpdcon", 
    "itmnormgpd", "itmweibullgpd", "itmgng", "mgammagpd", 
    "kden", "kdengpd", "kdengpdcon", "lognormgpd", "lognormgpdcon",
    "normgpd", "normgpdcon", "weibullgpd", "weibullgpdcon")

  if (!(modelname %in% allmodelnames))
    stop("invalid extreme value mixture model given")

  check.logic(upperfocus)
  check.logic(legend)
  check.prob(alpha, allownull = TRUE)
  if (!is.null(alpha)) 
    if ((alpha == 0) | (alpha == 1)) stop("alpha must be in (0, 1)")
  check.n(N)

  if (is.unsorted(modelfit$x, na.rm = TRUE)) {
    modelfit$x=sort(modelfit$x)
  } else {
    if (modelfit$x[1] > modelfit$x[length(modelfit$x)])
      modelfit$x = rev(modelfit$x)
  }
  
  nothresh = (modelname == "kden") | (modelname =="bckden") | (modelname =="dwm")
  twotail = (modelname == "gng") | (modelname == "gngcon") | (modelname == "gkg")| (modelname == "gkgcon")
  nophiu = (modelname == "itmnormgpd") | (modelname == "itmweibullgpd") | (modelname == "itmgng")
  
  if (nothresh) {
    upperfocus = FALSE
  }
  if (modelname == "gpd") {
    upperfocus = TRUE  
  }

  # makes following code easier to read
  if (twotail) {
    phiu = modelfit$phiur
    u = modelfit$ur
  } else if (nothresh) {
    phiu = -Inf
    u = Inf 
  } else if (nophiu) {
    phiu = -Inf
    u = ifelse(modelname == "itmgng", modelfit$ur, modelfit$u)
  } else {
    phiu = modelfit$phiu
    u = modelfit$u    
  }

  n = length(modelfit$x)
  emp.prob=ppoints(n)

  modelname = substring(strsplit(deparse(modelfit$call), "\\(")[[c(1, 1)]], 2)
  
  the.quant = switch(modelname,
    gpd = qgpd(emp.prob, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    bckden = qbckden(emp.prob, modelfit$kerncentres, modelfit$lambda, bw = NULL,
      modelfit$kernel, modelfit$bcmethod,
      modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
    bckdengpd = qbckdengpd(emp.prob, modelfit$x, modelfit$lambda, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel,
      modelfit$bcmethod, modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
    bckdengpdcon = qbckdengpdcon(emp.prob, modelfit$x, modelfit$lambda, modelfit$u,
      modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel,
      modelfit$bcmethod, modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
    betagpd = qbetagpd(emp.prob, modelfit$bshape1, modelfit$bshape2, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    betagpdcon = qbetagpdcon(emp.prob, modelfit$bshape1, modelfit$bshape2, modelfit$u,
      modelfit$xi, modelfit$phiu),
    dwm = qdwm(emp.prob, modelfit$wshape, modelfit$wscale, modelfit$cmu, modelfit$ctau,
      modelfit$sigmau, modelfit$xi),
    gammagpd = qgammagpd(emp.prob, modelfit$gshape, modelfit$gscale, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    gammagpdcon = qgammagpdcon(emp.prob, modelfit$gshape, modelfit$gscale, modelfit$u,
      modelfit$xi, modelfit$phiu),
    gkg = qgkg(emp.prob, modelfit$x, modelfit$lambda, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul,
      modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur, bw = NULL, modelfit$kernel),
    gkgcon = qgkgcon(emp.prob, modelfit$x, modelfit$lambda, 
      modelfit$ul, modelfit$xil, modelfit$phiul,
      modelfit$ur, modelfit$xir, modelfit$phiur, bw = NULL, modelfit$kernel),
    gng = qgng(emp.prob, modelfit$nmean, modelfit$nsd, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
      modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
    gngcon = qgngcon(emp.prob, modelfit$nmean, modelfit$nsd, 
      modelfit$ul, modelfit$xil, modelfit$phiul, modelfit$ur, modelfit$xir, modelfit$phiur),
    hpd = qhpd(emp.prob, modelfit$nmean, modelfit$nsd, modelfit$xi),
    hpdcon = qhpdcon(emp.prob, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$xi),
    itmnormgpd = qitmnormgpd(emp.prob, modelfit$nmean, modelfit$nsd, modelfit$epsilon, 
      modelfit$u, modelfit$sigmau, modelfit$xi),
    itmweibullgpd = qitmweibullgpd(emp.prob, modelfit$wshape, modelfit$wscale, modelfit$epsilon, 
      modelfit$u, modelfit$sigmau, modelfit$xi),
    itmgng = qitmgng(emp.prob, modelfit$nmean, modelfit$nsd, modelfit$epsilon, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$ur, modelfit$sigmaur, modelfit$xir),
    kden = qkden(emp.prob, modelfit$kerncentres, modelfit$lambda, bw = NULL, modelfit$kernel),
    kdengpd = qkdengpd(emp.prob, modelfit$x, modelfit$lambda, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel),
    kdengpdcon = qkdengpdcon(emp.prob, modelfit$x, modelfit$lambda, modelfit$u,
      modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel),
    mgammagpd = qmgammagpd(emp.prob, modelfit$mgshape, modelfit$mgscale, modelfit$mgweight,
      modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    lognormgpd = qlognormgpd(emp.prob, modelfit$lnmean, modelfit$lnsd, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    lognormgpdcon = qlognormgpd(emp.prob, modelfit$lnmean, modelfit$lnsd, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu), 
    normgpd = qnormgpd(emp.prob, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau,
      modelfit$xi, modelfit$phiu),
    normgpdcon = qnormgpdcon(emp.prob, modelfit$nmean, modelfit$nsd, modelfit$u,
      modelfit$xi, modelfit$phiu),
    weibullgpd = qweibullgpd(emp.prob, modelfit$wshape, modelfit$wscale, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    weibullgpdcon = qweibullgpdcon(emp.prob, modelfit$wshape, modelfit$wscale, modelfit$u,
      modelfit$xi, modelfit$phiu))

  axislims = range(c(modelfit$x, the.quant), na.rm = TRUE)

  if (!is.null(alpha)) {
    # obtain CI by Monte Carlo simulation (ignores estimation uncertainty)
    simulate.new.quantiles <- function(i, modelfit, modelname) {
      simdata = switch(modelname,
        gpd = rgpd(n, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        bckden = rbckden(n, modelfit$kerncentres, modelfit$lambda, bw = NULL, modelfit$kernel,
          modelfit$bcmethod, modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
        bckdengpd = rbckdengpd(n, modelfit$x, modelfit$lambda, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel,
          modelfit$bcmethod, modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
        bckdengpdcon = rbckdengpdcon(n, modelfit$x, modelfit$lambda, modelfit$u,
          modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel,
          modelfit$bcmethod, modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
        betagpd = rbetagpd(n, modelfit$bshape1, modelfit$bshape2, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        betagpdcon = rbetagpdcon(n, modelfit$bshape1, modelfit$bshape2, modelfit$u,
          modelfit$xi, modelfit$phiu),
        dwm = rdwm(n, modelfit$wshape, modelfit$wscale, modelfit$cmu, modelfit$ctau,
          modelfit$sigmau, modelfit$xi),
        gammagpd = rgammagpd(n, modelfit$gshape, modelfit$gscale, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        gammagpdcon = rgammagpdcon(n, modelfit$gshape, modelfit$gscale, modelfit$u,
          modelfit$xi, modelfit$phiu),
        gkg = rgkg(n, modelfit$x, modelfit$lambda, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul,
          modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur, bw = NULL, modelfit$kernel),
        gkgcon = rgkgcon(n, modelfit$x, modelfit$lambda, 
          modelfit$ul, modelfit$xil, modelfit$phiul,
          modelfit$ur, modelfit$xir, modelfit$phiur, bw = NULL, modelfit$kernel),
        gng = rgng(n, modelfit$nmean, modelfit$nsd, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
          modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
        gngcon = rgngcon(n, modelfit$nmean, modelfit$nsd, 
          modelfit$ul, modelfit$xil, modelfit$phiul,
          modelfit$ur, modelfit$xir, modelfit$phiur),
        hpd = rhpd(n, modelfit$nmean, modelfit$nsd, modelfit$xi),
        hpdcon = rhpdcon(n, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$xi),
        itmnormgpd = ritmnormgpd(n, modelfit$nmean, modelfit$nsd, modelfit$epsilon, 
          modelfit$u, modelfit$sigmau, modelfit$xi),
        itmweibullgpd = ritmweibullgpd(n, modelfit$wshape, modelfit$wscale, modelfit$epsilon, 
          modelfit$u, modelfit$sigmau, modelfit$xi),
        itmgng = ritmgng(n, modelfit$nmean, modelfit$nsd, modelfit$epsilon, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$ur, modelfit$sigmaur, modelfit$xir),
        kden = rkden(n, modelfit$kerncentres, modelfit$lambda, bw = NULL, modelfit$kernel),
        kdengpd = rkdengpd(n, modelfit$x, modelfit$lambda, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel),
        kdengpdcon = rkdengpdcon(n, modelfit$x, modelfit$lambda, modelfit$u,
          modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel),
        lognormgpd = rlognormgpd(n, modelfit$lnmean, modelfit$lnsd, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        lognormgpdcon = rlognormgpd(n, modelfit$lnmean, modelfit$lnsd, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        mgammagpd = rmgammagpd(n, modelfit$mgshape, modelfit$mgscale, modelfit$mgweight,
          modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        normgpd = rnormgpd(n, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau,
          modelfit$xi, modelfit$phiu),
        normgpdcon = rnormgpdcon(n, modelfit$nmean, modelfit$nsd, modelfit$u,
          modelfit$xi, modelfit$phiu),
        weibullgpd = rweibullgpd(n, modelfit$wshape, modelfit$wscale, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        weibullgpdcon = rweibullgpdcon(n, modelfit$wshape, modelfit$wscale, modelfit$u,
          modelfit$xi, modelfit$phiu))
      sort(simdata, na.last = FALSE)
    }
    sim.q = sapply(1:N, FUN = simulate.new.quantiles, modelfit = modelfit, modelname = modelname)
    ci.q = apply(sim.q, FUN = quantile, MARGIN = 1, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
    # simulation for gpd allow for uncertainty in phiu, so need to ignore CI below threshold
    if (modelname == "gpd") {
      ci.q[1, modelfit$x < u] = NA
      ci.q[2, modelfit$x < u] = NA
    }
    axislims[2] = max(c(ci.q, axislims[2]), na.rm = TRUE)
    axislims[1] = min(c(ci.q, axislims[1]), na.rm = TRUE)
  }

  if (upperfocus) {
    axislims = c(u, axislims[2])    
  }
  
  if (is.infinite(axislims[2])) axislims[2] = max(modelfit$x)

  plot(the.quant, modelfit$x, pch = 'x',
    xlab = "Theoretical Quantiles", ylab = "Empirical Quantiles", 
    xlim = axislims, ylim = axislims, ...)
  abline(c(0, 1), col = "red")

  # add threshold indicator lines
  if ((twotail) | (modelname == "itmgng")) {
    abline(h = modelfit$ul, v = modelfit$ul, lty = 3)
    abline(h = modelfit$ur, v = modelfit$ur, lty = 3)
  } else {
    abline(h = u, v = u, lty = 3)
  }
  
  if (!is.null(alpha)) {
    lines(the.quant, ci.q[1, ], lty=2)
    lines(the.quant, ci.q[2, ], lty=2)
    if (legend) {
      if (nothresh) {
        legend("bottomright", c("Data", "Line of Equality",
          paste("Simulated Pointwise ", formatC(100 * (1 - alpha), format = "f", digits=1), "% CI", sep="")),
          pch = c(1, rep(-1, 2)), lty = c(0, 1, 2), bg = "white", col = c("black", "red", "black"))        
      } else {
        legend("bottomright", c("Data", "Line of Equality",
          paste("Simulated Pointwise ", formatC(100 * (1 - alpha), format = "f", digits=1), "% CI", sep=""),
          ifelse(((twotail) | (modelname == "itmgng")) & !upperfocus, "Thresholds", "Threshold")),
          pch = c(1, rep(-1, 3)), lty = c(0, 1, 2, 3), bg = "white", col = c("black", "red", rep("black", 2)))
      }
    }
  } else {
    if (legend) {
      if (nothresh) {
        legend("bottomright", c("Data", "Line of Equality"),
          pch = c(1, -1), lty = c(0, 1), bg = "white", col = c("black", "red", "black"))  
      } else {
        legend("bottomright", c("Data", "Line of Equality",
          ifelse(((twotail) | (modelname == "itmgng")) & !upperfocus, "Thresholds" ,"Threshold")),
          pch = c(1, rep(-1, 2)), lty = c(0, 1, 3), bg = "white", col = c("black", "red", "black"))
      }
    }
  }
}

#' @export
#' @aliases evmix.diag rlplot qplot pplot densplot
#' @rdname evmix.diag

pplot <- function(modelfit, upperfocus = TRUE, alpha = 0.05, N = 1000,
  legend = TRUE, ...) {
  
  # Check properties of inputs
  if (missing(modelfit))
    stop("modelfit from extreme value mixture model must be given")
   
  if (!exists("call", modelfit))
    stop("modelfit from extreme value mixture model must be given")

  modelname = substring(strsplit(deparse(modelfit$call),"\\(")[[c(1, 1)]], 2)

  allmodelnames = c("gpd", "betagpd", "betagpdcon", 
    "bckden", "bckdengpd", "bckdengpdcon", "dwm", "gammagpd", "gammagpdcon",
    "gkg", "gkgcon", "gng", "gngcon", "hpd", "hpdcon", 
    "itmnormgpd", "itmweibullgpd", "itmgng", "mgammagpd", 
    "kden", "kdengpd", "kdengpdcon", "lognormgpd", "lognormgpdcon",
    "normgpd", "normgpdcon", "weibullgpd", "weibullgpdcon")

  if (!(modelname %in% allmodelnames))
    stop("invalid extreme value mixture model given")

  check.logic(upperfocus)
  check.logic(legend)
  check.prob(alpha, allownull = TRUE)
  if (!is.null(alpha)) 
    if ((alpha == 0) | (alpha == 1)) stop("alpha must be in (0, 1)")
  check.n(N)

  if (is.unsorted(modelfit$x, na.rm = TRUE)) {
    modelfit$x=sort(modelfit$x)
  } else {
    if (modelfit$x[1] > modelfit$x[length(modelfit$x)])
      modelfit$x = rev(modelfit$x)
  }
  
  nothresh = (modelname == "kden") | (modelname =="bckden") | (modelname =="dwm")
  twotail = (modelname == "gng") | (modelname == "gngcon") | (modelname == "gkg") | (modelname == "gkgcon")
  nophiu = (modelname == "itmnormgpd") | (modelname == "itmweibullgpd") | (modelname == "itmgng")
  
  if (nothresh) {
    upperfocus = FALSE
  }
  if (modelname == "gpd") {
    upperfocus = TRUE  
  }

  # makes following code easier to read
  if (twotail) {
    phiu = modelfit$phiur
    u = modelfit$ur
  } else if (nothresh) {
    phiu = -Inf
    u = Inf
  } else if (nophiu) {
    phiu = -Inf
    u = ifelse(modelname == "itmgng", modelfit$ur, modelfit$u)
  } else {
    phiu = modelfit$phiu
    u = modelfit$u    
  }

  n = length(modelfit$x)
  emp.prob = ppoints(n)

  the.prob = switch(modelname,
    gpd = pgpd(modelfit$x, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    bckden = pbckden(modelfit$x, modelfit$kerncentres, modelfit$lambda, bw = NULL, modelfit$kernel,
      modelfit$bcmethod, modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
    bckdengpd = pbckdengpd(modelfit$x, modelfit$x, modelfit$lambda, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel, 
      modelfit$bcmethod, modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
    bckdengpdcon = pbckdengpdcon(modelfit$x, modelfit$x, modelfit$lambda, modelfit$u,
      modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel, 
      modelfit$bcmethod, modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
    betagpd = pbetagpd(modelfit$x, modelfit$bshape1, modelfit$bshape2, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    betagpdcon = pbetagpdcon(modelfit$x, modelfit$bshape1, modelfit$bshape2, modelfit$u,
      modelfit$xi, modelfit$phiu),
    dwm = pdwm(modelfit$x, modelfit$wshape, modelfit$wscale, modelfit$cmu, modelfit$ctau,
      modelfit$sigmau, modelfit$xi),
    gammagpd = pgammagpd(modelfit$x, modelfit$gshape, modelfit$gscale, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    gammagpdcon = pgammagpdcon(modelfit$x, modelfit$gshape, modelfit$gscale, modelfit$u,
      modelfit$xi, modelfit$phiu),
    gkg = pgkg(modelfit$x, modelfit$x, modelfit$lambda, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul,
      modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur, bw = NULL, modelfit$kernel),
    gkgcon = pgkgcon(modelfit$x, modelfit$x, modelfit$lambda, 
      modelfit$ul, modelfit$xil, modelfit$phiul,
      modelfit$ur, modelfit$xir, modelfit$phiur, bw = NULL, modelfit$kernel),
    gng = pgng(modelfit$x, modelfit$nmean, modelfit$nsd,
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
      modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
    gngcon = pgngcon(modelfit$x, modelfit$nmean, modelfit$nsd, 
      modelfit$ul, modelfit$xil, modelfit$phiul, modelfit$ur, modelfit$xir, modelfit$phiur),
    hpd = phpd(modelfit$x, modelfit$nmean, modelfit$nsd, modelfit$xi),
    hpdcon = phpdcon(modelfit$x, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$xi),
    itmnormgpd = pitmnormgpd(modelfit$x, modelfit$nmean, modelfit$nsd, modelfit$epsilon, 
      modelfit$u, modelfit$sigmau, modelfit$xi),
    itmweibullgpd = pitmweibullgpd(modelfit$x, modelfit$wshape, modelfit$wscale, modelfit$epsilon, 
      modelfit$u, modelfit$sigmau, modelfit$xi),
    itmgng = pitmgng(modelfit$x, modelfit$nmean, modelfit$nsd, modelfit$epsilon, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$ur, modelfit$sigmaur, modelfit$xir),
    kden = pkden(modelfit$x, modelfit$kerncentres, modelfit$lambda, bw = NULL, modelfit$kernel),
    kdengpd = pkdengpd(modelfit$x, modelfit$x, modelfit$lambda, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel),
    kdengpdcon = pkdengpdcon(modelfit$x, modelfit$x, modelfit$lambda, modelfit$u,
      modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel),
    lognormgpd = plognormgpd(modelfit$x, modelfit$lnmean, modelfit$lnsd, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    lognormgpdcon = plognormgpd(modelfit$x, modelfit$lnmean, modelfit$lnsd, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    mgammagpd = pmgammagpd(modelfit$x, modelfit$mgshape, modelfit$mgscale, modelfit$mgweight,
      modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    normgpd = pnormgpd(modelfit$x, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau,
      modelfit$xi, modelfit$phiu),
    normgpdcon = pnormgpdcon(modelfit$x, modelfit$nmean, modelfit$nsd, modelfit$u,
      modelfit$xi, modelfit$phiu),
    weibullgpd = pweibullgpd(modelfit$x, modelfit$wshape, modelfit$wscale, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    weibullgpdcon = pweibullgpdcon(modelfit$x, modelfit$wshape, modelfit$wscale, modelfit$u,
      modelfit$xi, modelfit$phiu))

  # If focus on upper tail then must display at least upper 10%
  # even if tail fraction is smaller, otherwise may not look nice
  if (upperfocus) { 
    axislims = c(1 - ifelse(phiu < 0.1, 0.1, phiu), 1)
  } else {
    axislims = c(0, 1)
  }

  if (!is.null(alpha)) {
    # obtain CI by Monte Carlo simulation (ignores estimation uncertainty)
    simulate.new.probabilities <- function(i, modelfit, modelname) {
      simdata = switch(modelname,
        gpd = rgpd(n, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        bckden = rbckden(n, modelfit$kerncentres, modelfit$lambda, bw = NULL, modelfit$kernel,
          modelfit$bcmethod, modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
        bckdengpd = rbckdengpd(n, modelfit$x, modelfit$lambda, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel,
          modelfit$bcmethod, modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
        bckdengpdcon = rbckdengpdcon(n, modelfit$x, modelfit$lambda, modelfit$u,
          modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel,
          modelfit$bcmethod, modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
        betagpd = rbetagpd(n, modelfit$bshape1, modelfit$bshape2, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        betagpdcon = rbetagpdcon(n, modelfit$bshape1, modelfit$bshape2, modelfit$u,
          modelfit$xi, modelfit$phiu),
        dwm = rdwm(n, modelfit$wshape, modelfit$wscale, modelfit$cmu, modelfit$ctau,
          modelfit$sigmau, modelfit$xi),
        gammagpd = rgammagpd(n, modelfit$gshape, modelfit$gscale, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        gammagpdcon = rgammagpdcon(n, modelfit$gshape, modelfit$gscale, modelfit$u,
          modelfit$xi, modelfit$phiu),
        gkg = rgkg(n, modelfit$x, modelfit$lambda, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul,
          modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur, bw = NULL, modelfit$kernel),
        gkgcon = rgkgcon(n, modelfit$x, modelfit$lambda, 
          modelfit$ul, modelfit$xil, modelfit$phiul,
          modelfit$ur, modelfit$xir, modelfit$phiur, bw = NULL, modelfit$kernel),
        gng = rgng(n, modelfit$nmean, modelfit$nsd, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
          modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
        gngcon = rgngcon(n, modelfit$nmean, modelfit$nsd, 
          modelfit$ul, modelfit$xil, modelfit$phiul,
          modelfit$ur, modelfit$xir, modelfit$phiur),
        hpd = rhpd(n, modelfit$nmean, modelfit$nsd, modelfit$xi),
        hpdcon = rhpdcon(n, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$xi),
        itmnormgpd = ritmnormgpd(n, modelfit$nmean, modelfit$nsd, modelfit$epsilon, 
          modelfit$u, modelfit$sigmau, modelfit$xi),
        itmweibullgpd = ritmweibullgpd(n, modelfit$wshape, modelfit$wscale, modelfit$epsilon, 
          modelfit$u, modelfit$sigmau, modelfit$xi),
        itmgng = ritmgng(n, modelfit$nmean, modelfit$nsd, modelfit$epsilon, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$ur, modelfit$sigmaur, modelfit$xir),
        kden = rkden(n, modelfit$kerncentres, modelfit$lambda, bw = NULL, modelfit$kernel),
        kdengpd = rkdengpd(n, modelfit$x, modelfit$lambda, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel),
        kdengpdcon = rkdengpdcon(n, modelfit$x, modelfit$lambda, modelfit$u,
          modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel),
        lognormgpd = rlognormgpd(n, modelfit$lnmean, modelfit$lnsd, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        lognormgpdcon = rlognormgpd(n, modelfit$lnmean, modelfit$lnsd, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        mgammagpd = rmgammagpd(n, modelfit$mgshape, modelfit$mgscale, modelfit$mgweight,
          modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        normgpd = rnormgpd(n, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau,
          modelfit$xi, modelfit$phiu),
        normgpdcon = rnormgpdcon(n, modelfit$nmean, modelfit$nsd, modelfit$u,
          modelfit$xi, modelfit$phiu),
        weibullgpd = rweibullgpd(n, modelfit$wshape, modelfit$wscale, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        weibullgpdcon = rweibullgpdcon(n, modelfit$wshape, modelfit$wscale, modelfit$u,
          modelfit$xi, modelfit$phiu))
      
      simdata = sort(simdata, na.last = FALSE)
      simprob = switch(modelname,
        gpd = pgpd(simdata, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        bckden = pbckden(simdata, modelfit$kerncentres, modelfit$lambda, bw = NULL, modelfit$kernel,
          modelfit$bcmethod, modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
        bckdengpd = pbckdengpd(simdata, modelfit$x, modelfit$lambda, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel,
          modelfit$bcmethod, modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),  
        bckdengpdcon = pbckdengpdcon(simdata, modelfit$x, modelfit$lambda, modelfit$u,
          modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel, modelfit$bcmethod,
          modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),  
        betagpd = pbetagpd(simdata, modelfit$bshape1, modelfit$bshape2, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        betagpdcon = pbetagpdcon(simdata, modelfit$bshape1, modelfit$bshape2, modelfit$u,
          modelfit$xi, modelfit$phiu),
        dwm = pdwm(simdata, modelfit$wshape, modelfit$wscale, modelfit$cmu, modelfit$ctau,
          modelfit$sigmau, modelfit$xi),
        gammagpd = pgammagpd(simdata, modelfit$gshape, modelfit$gscale, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        gammagpdcon = pgammagpdcon(simdata, modelfit$gshape, modelfit$gscale, modelfit$u,
          modelfit$xi, modelfit$phiu),
        gkg = pgkg(simdata, modelfit$x, modelfit$lambda, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul,
          modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur, bw = NULL, modelfit$kernel),
        gkgcon = pgkgcon(simdata, modelfit$x, modelfit$lambda, 
          modelfit$ul, modelfit$xil, modelfit$phiul,
          modelfit$ur, modelfit$xir, modelfit$phiur, bw = NULL, modelfit$kernel),
        gng = pgng(simdata, modelfit$nmean, modelfit$nsd, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
          modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
        gngcon = pgngcon(simdata, modelfit$nmean, modelfit$nsd, 
          modelfit$ul, modelfit$xil, modelfit$phiul, modelfit$ur, modelfit$xir, modelfit$phiur),
        hpd = phpd(simdata, modelfit$nmean, modelfit$nsd, modelfit$xi),
        hpdcon = phpdcon(simdata, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$xi),
        itmnormgpd = pitmnormgpd(simdata, modelfit$nmean, modelfit$nsd, modelfit$epsilon, 
          modelfit$u, modelfit$sigmau, modelfit$xi),
        itmweibullgpd = pitmweibullgpd(simdata, modelfit$wshape, modelfit$wscale, modelfit$epsilon, 
          modelfit$u, modelfit$sigmau, modelfit$xi),
        itmgng = pitmgng(simdata, modelfit$nmean, modelfit$nsd, modelfit$epsilon, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$ur, modelfit$sigmaur, modelfit$xir),
        kden = pkden(simdata, modelfit$kerncentres, modelfit$lambda, bw = NULL, modelfit$kernel),
        kdengpd = pkdengpd(simdata, modelfit$x, modelfit$lambda, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel),
        kdengpdcon = pkdengpdcon(simdata, modelfit$x, modelfit$lambda, modelfit$u,
          modelfit$xi, modelfit$phiu, bw = NULL, modelfit$kernel),
        lognormgpd = plognormgpd(simdata, modelfit$lnmean, modelfit$lnsd, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        lognormgpdcon = plognormgpd(simdata, modelfit$lnmean, modelfit$lnsd, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        mgammagpd = pmgammagpd(simdata, modelfit$mgshape, modelfit$mgscale, modelfit$mgweight,
          modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        normgpd = pnormgpd(simdata, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau,
          modelfit$xi, modelfit$phiu),
        normgpdcon = pnormgpdcon(simdata, modelfit$nmean, modelfit$nsd, modelfit$u,
          modelfit$xi, modelfit$phiu),
        weibullgpd = pweibullgpd(simdata, modelfit$wshape, modelfit$wscale, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        weibullgpdcon = pweibullgpdcon(simdata, modelfit$wshape, modelfit$wscale, modelfit$u,
          modelfit$xi, modelfit$phiu))
      simprob      
    }
    sim.p = sapply(1:N, FUN = simulate.new.probabilities, modelfit = modelfit, modelname = modelname)
    ci.p = apply(sim.p, FUN = quantile, MARGIN = 1, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  }

  plot(the.prob, emp.prob, pch = 'x', cex = 0.8,
    xlab = "Theoretical Probability", ylab = "Empirical Probability",
    xlim = axislims, ylim = axislims, ...)
  abline(c(0, 1), col = "red")

  # add tail fraction indicator lines
  if (twotail) {
    abline(h = modelfit$phiul, v = modelfit$phiul, lty = 3)
    abline(h = 1 - modelfit$phiur, v = 1 - modelfit$phiur, lty = 3)
  } else {
    abline(h = 1 - phiu, v = 1 - phiu, lty = 3)
  }
  
  if (!is.null(alpha)) {
    lines(the.prob, ci.p[1, ], lty=2)
    lines(the.prob, ci.p[2, ], lty=2)
    if (legend) {
      if (nophiu) {
        legend("bottomright", c("Data", "Line of Equality",
          paste("Simulated Pointwise ", formatC(100 * (1 - alpha), format = "f", digits=1), "% CI", sep="")),
          pch = c(1, rep(-1, 2)), lty = c(0, 1, 2), bg = "white", col = c("black", "red", "black"))
      } else {
        legend("bottomright", c("Data", "Line of Equality",
          paste("Simulated Pointwise ", formatC(100 * (1 - alpha), format = "f", digits=1), "% CI", sep=""),
          ifelse(twotail & !upperfocus, "Tail Fractions", "Tail Fraction")),
          pch = c(1, rep(-1, 3)), lty = c(0, 1, 2, 3), bg = "white", col = c("black", "red", rep("black", 2)))        
      }
    }
  } else {
    if (legend) {
      if (nophiu) {
        legend("bottomright", c("Data", "Line of Equality"),
          pch = c(1, -1), lty = c(0, 1), bg = "white", col = c("black", "red"))
      } else {
        legend("bottomright", c("Data", "Line of Equality",
          ifelse(twotail & !upperfocus, "Tail Fractions", "Tail Fraction")),
          pch = c(1, rep(-1, 2)), lty = c(0, 1, 3), bg = "white", col = c("black", "red", "black"))
      }
    }
  }
}

#' @export
#' @aliases evmix.diag rlplot qplot pplot densplot
#' @rdname evmix.diag

densplot <- function(modelfit, upperfocus = TRUE, legend = TRUE, ...) {
    
  # Check properties of inputs
  if (missing(modelfit))
    stop("modelfit from extreme value mixture model must be given")
   
  if (!exists("call", modelfit))
    stop("modelfit from extreme value mixture model must be given")

  modelname = substring(strsplit(deparse(modelfit$call), "\\(")[[c(1, 1)]], 2)

  allmodelnames = c("gpd", "betagpd", "betagpdcon", 
    "bckden", "bckdengpd", "bckdengpdcon", "dwm", "gammagpd", "gammagpdcon",
    "gkg", "gkgcon", "gng", "gngcon", "hpd", "hpdcon", 
    "itmnormgpd", "itmweibullgpd", "itmgng", "mgammagpd", 
    "kden", "kdengpd", "kdengpdcon", "lognormgpd", "lognormgpdcon",
    "normgpd", "normgpdcon", "weibullgpd", "weibullgpdcon")

  if (!(modelname %in% allmodelnames))
    stop("invalid extreme value mixture model given")

  check.logic(upperfocus)
  check.logic(legend)
  
  if (is.unsorted(modelfit$x, na.rm = TRUE)) {
    modelfit$x=sort(modelfit$x)
  } else {
    if (modelfit$x[1] > modelfit$x[length(modelfit$x)])
      modelfit$x = rev(modelfit$x)
  }

  nothresh = (modelname == "kden") | (modelname =="bckden") | (modelname =="dwm")
  twotail = (modelname == "gng") | (modelname == "gngcon") | (modelname == "gkg") | (modelname == "gkgcon")
  
  nophiu = (modelname == "itmnormgpd") | (modelname == "itmweibullgpd") | (modelname == "itmgng")
  
  if (nothresh) {
    upperfocus = FALSE
  }
  if (modelname == "gpd") {
    upperfocus = TRUE  
  }

  # makes following code easier to read
  if (twotail) {
    phiu = modelfit$phiur
    u = modelfit$ur
  } else if (nothresh) {
    phiu = -Inf
    u = Inf 
  } else if (nophiu) {
    phiu = -Inf
    u = ifelse(modelname == "itmgng", modelfit$ur, modelfit$u)
  } else {
    phiu = modelfit$phiu
    u = modelfit$u    
  }
  
  xx = seq(min(modelfit$x, na.rm = TRUE), max(modelfit$x, na.rm = TRUE), sd(modelfit$x, na.rm = TRUE)/100)
  
  the.dens = switch(modelname,
    gpd = dgpd(xx, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    bckden = dbckden(xx, modelfit$kerncentres, modelfit$lambda, bw=NULL, modelfit$kernel,
      modelfit$bcmethod, modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
    bckdengpd = dbckdengpd(xx, modelfit$x, modelfit$lambda, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu, bw=NULL, modelfit$kernel,
      modelfit$bcmethod, modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
    bckdengpdcon = dbckdengpdcon(xx, modelfit$x, modelfit$lambda, modelfit$u,
      modelfit$xi, modelfit$phiu, bw=NULL, modelfit$kernel, modelfit$bcmethod,
      modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
    betagpd = dbetagpd(xx, modelfit$bshape1, modelfit$bshape2, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    betagpdcon = dbetagpdcon(xx, modelfit$bshape1, modelfit$bshape2, modelfit$u,
      modelfit$xi, modelfit$phiu),
    dwm = ddwm(xx, modelfit$wshape, modelfit$wscale, modelfit$cmu, modelfit$ctau,
      modelfit$sigmau, modelfit$xi),
    gammagpd = dgammagpd(xx, modelfit$gshape, modelfit$gscale, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    gammagpdcon = dgammagpdcon(xx, modelfit$gshape, modelfit$gscale, modelfit$u,
      modelfit$xi, modelfit$phiu),
    gkg = dgkg(xx, modelfit$x, modelfit$lambda, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul,
      modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur, bw=NULL, modelfit$kernel),
    gkgcon = dgkgcon(xx, modelfit$x, modelfit$lambda, 
      modelfit$ul, modelfit$xil, modelfit$phiul,
      modelfit$ur, modelfit$xir, modelfit$phiur, bw=NULL, modelfit$kernel),
    gng = dgng(xx, modelfit$nmean, modelfit$nsd, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
      modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
    gngcon = dgngcon(xx, modelfit$nmean, modelfit$nsd, 
      modelfit$ul, modelfit$xil, modelfit$phiul, modelfit$ur, modelfit$xir, modelfit$phiur),
    hpd = dhpd(xx, modelfit$nmean, modelfit$nsd, modelfit$xi),
    hpdcon = dhpdcon(xx, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$xi),
    itmnormgpd = ditmnormgpd(xx, modelfit$nmean, modelfit$nsd, modelfit$epsilon, 
      modelfit$u, modelfit$sigmau, modelfit$xi),
    itmweibullgpd = ditmweibullgpd(xx, modelfit$wshape, modelfit$wscale, modelfit$epsilon, 
      modelfit$u, modelfit$sigmau, modelfit$xi),
    itmgng = ditmgng(xx, modelfit$nmean, modelfit$nsd, modelfit$epsilon, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$ur, modelfit$sigmaur, modelfit$xir),
    kden = dkden(xx, modelfit$kerncentres, modelfit$lambda, bw=NULL, modelfit$kernel),
    kdengpd = dkdengpd(xx, modelfit$x, modelfit$lambda, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu, bw=NULL, modelfit$kernel),
    kdengpdcon = dkdengpdcon(xx, modelfit$x, modelfit$lambda, modelfit$u,
      modelfit$xi, modelfit$phiu, bw=NULL, modelfit$kernel),
    lognormgpd = dlognormgpd(xx, modelfit$lnmean, modelfit$lnsd, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    lognormgpdcon = dlognormgpd(xx, modelfit$lnmean, modelfit$lnsd, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    mgammagpd = dmgammagpd(xx, modelfit$mgshape, modelfit$mgscale, modelfit$mgweight,
      modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    normgpd = dnormgpd(xx, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau,
      modelfit$xi, modelfit$phiu),
    normgpdcon = dnormgpdcon(xx, modelfit$nmean, modelfit$nsd, modelfit$u,
      modelfit$xi, modelfit$phiu),
    weibullgpd = dweibullgpd(xx, modelfit$wshape, modelfit$wscale, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    weibullgpdcon = dweibullgpdcon(xx, modelfit$wshape, modelfit$wscale, modelfit$u,
      modelfit$xi, modelfit$phiu))

  kde.est = density(modelfit$x, from = min(modelfit$x, na.rm = TRUE),
    to = max(modelfit$x, na.rm = TRUE), na.rm = TRUE)
  
  # In GPD case NA are treated as below threshold and retained from fitting
  # so need to scale KDE
  if ((modelname == "gpd") & any(is.na(modelfit$x))) {
    kde.est$y = kde.est$y * modelfit$phiu
  }
  
  if (upperfocus) {
    xlims = c(u, max(modelfit$x, na.rm = TRUE))    
    ylims = c(0, 1.4 * max(c(the.dens[xx > u], kde.est$y[kde.est$x > u])))    
  } else {
    xlims = range(modelfit$x, na.rm = TRUE)
    ylims = c(0, 1.1 * max(c(the.dens, kde.est$y)))    
  }
  
  if ((modelname == "gpd") & any(is.na(modelfit$x))) {
    hist(ifelse(is.na(modelfit$x), u - (max(modelfit$x, na.rm = TRUE) - u), modelfit$x),
      breaks = "Scott", freq = FALSE, xlab = "Sample Data", main="", xlim = xlims, ylim = ylims)
  } else {
    hist(modelfit$x, breaks = "Scott", freq = FALSE, 
      xlab = "Sample Data", main="", xlim = xlims, ylim = ylims)    
  }
  rug(modelfit$x, quiet = TRUE)
  box()
    
  lines(xx, the.dens, lwd = 2)
  lines(kde.est, lty = 2, lwd = 1.5, col = "green")

  # add threshold indicator lines
  if ((twotail) | (modelname == "itmgng")) {
    abline(v = c(modelfit$ul, modelfit$ur), lty = 3)
  } else {
    abline(v = u, lty = 3)
  }
  
  if (legend) {
    if (nothresh){
      legend("topright", c("Fitted Density", "KDE"),
        lty = c(1, 2), lwd = c(2, 1.5), bg = "white", col = c("black", "green"))
    } else {
      legend("topright", c("Fitted Density", 
        ifelse((twotail) | (modelname == "itmgng"), "Thresholds", "Threshold"), "KDE"),
        lty = c(1, 3, 2), lwd = c(2, 1, 1.5), bg = "white", col = c(rep("black", 2), "green"))      
    }
  }
}
