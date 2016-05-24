#' @name dwm
#' 
#' @title Dynamically Weighted Mixture Model
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the dynamically weighted mixture model. The
#'   parameters are the Weibull shape \code{wshape} and scale \code{wscale},
#'   Cauchy location \code{cmu}, Cauchy scale \code{ctau}, GPD scale
#'   \code{sigmau}, shape \code{xi} and initial value for the quantile
#'   \code{qinit}.
#'
#' @inheritParams weibullgpd
#' @param cmu     Cauchy location
#' @param ctau    Cauchy scale
#' @param qinit   scalar or vector of initial values for the quantile estimate
#' @inheritParams gpd
#' 
#' @details The dynamic weighted mixture model combines a Weibull for the bulk 
#'   model with GPD for the tail model. However, unlike all the other mixture 
#'   models the GPD is defined over the entire range of support rather than as a
#'   conditional model above some threshold. A transition function is used to 
#'   apply weights to transition between the bulk and GPD for the upper tail, 
#'   thus providing the dynamically weighted mixture. They use a Cauchy 
#'   cumulative distribution function for the transition function.
#'   
#'   The density function is then a dynamically weighted mixture given by: 
#'   \deqn{f(x) = {[1 - p(x)] h(x) + p(x) g(x)}/r} where \eqn{h(x)} and
#'   \eqn{g(x)} are the Weibull and unscaled GPD density functions respectively
#'   (i.e. \code{dweibull(x, wshape, wscale)} and \code{dgpd(x, u, sigmau,
#'   xi)}). The Cauchy cumulative distribution function used to provide the
#'   transition is defined by \eqn{p(x)} (i.e. \code{pcauchy(x, cmu, ctau}. The
#'   normalisation constant \eqn{r} ensures a proper density.
#'   
#'   The quantile function is not available in closed form, so has to be solved 
#'   numerically. The argument \code{qinit} is the initial quantile estimate
#'   which is used for numerical optimisation and should be set to a reasonable
#'   guess. When the \code{qinit} is \code{NULL}, the initial quantile value is
#'   given by the midpoint between the Weibull and GPD quantiles. As with the
#'   other inputs \code{qinit} is also vectorised, but \code{R} does not permit
#'   vectors combining \code{NULL} and numeric entries.
#' 
#' @return \code{\link[evmix:dwm]{ddwm}} gives the density, 
#' \code{\link[evmix:dwm]{pdwm}} gives the cumulative distribution function,
#' \code{\link[evmix:dwm]{qdwm}} gives the quantile function and 
#' \code{\link[evmix:dwm]{rdwm}} gives a random sample.
#' 
#' @note All inputs are vectorised except \code{log} and \code{lower.tail}.
#' The main inputs (\code{x}, \code{p} or \code{q}) and parameters must be either
#' a scalar or a vector. If vectors are provided they must all be of the same length,
#' and the function will be evaluated for each element of vector. In the case of 
#' \code{\link[evmix:dwm]{rdwm}} any input vector must be of length \code{n}.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:dwm]{rdwm}} is 1.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x},
#' \code{p} and \code{q} are passed through as is and infinite values are set to
#' \code{NA}. None of these are not permitted for the parameters.
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' \url{http://en.wikipedia.org/wiki/Weibull_distribution}
#' 
#' \url{http://en.wikipedia.org/wiki/Cauchy_distribution}
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Frigessi, A., Haug, O. and Rue, H. (2002). A dynamic mixture model for unsupervised tail
#' estimation without threshold selection. Extremes 5 (3), 219-235
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:gpd]{gpd}}, \code{\link[stats:Cauchy]{dcauchy}}
#'    and \code{\link[stats:Weibull]{dweibull}}
#' @aliases dwm ddwm pdwm qdwm rdwm
#' @family  ldwm fdwm
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(2, 2))
#' 
#' xx = seq(0.001, 5, 0.01)
#' f = ddwm(xx, wshape = 2, wscale = 1/gamma(1.5), cmu = 1, ctau = 1, sigmau = 1, xi = 0.5)
#' plot(xx, f, ylim = c(0, 1), xlim = c(0, 5), type = 'l', lwd = 2, 
#'   ylab = "density", main = "Plot example in Frigessi et al. (2002)")
#' lines(xx, dgpd(xx, sigmau = 1, xi = 0.5), col = "red", lty = 2, lwd = 2)
#' lines(xx, dweibull(xx, shape = 2, scale = 1/gamma(1.5)), col = "blue", lty = 2, lwd = 2)
#' legend('topright', c('DWM', 'Weibull', 'GPD'),
#'       col = c("black", "blue", "red"), lty = c(1, 2, 2), lwd = 2)
#' 
#' # three tail behaviours
#' plot(xx, pdwm(xx, xi = 0), type = "l")
#' lines(xx, pdwm(xx, xi = 0.3), col = "red")
#' lines(xx, pdwm(xx, xi = -0.3), col = "blue")
#' legend("bottomright", paste("xi =",c(0, 0.3, -0.3)), col=c("black", "red", "blue"), lty = 1)
#' 
#' x = rdwm(10000, wshape = 2, wscale = 1/gamma(1.5), cmu = 1, ctau = 1, sigmau = 1, xi = 0.1)
#' xx = seq(0, 15, 0.01)
#' hist(x, freq = FALSE, breaks = 100)
#' lines(xx, ddwm(xx, wshape = 2, wscale = 1/gamma(1.5), cmu = 1, ctau = 1, sigmau = 1, xi = 0.1),
#'   lwd = 2, col = 'black')
#'   
#' plot(xx, pdwm(xx, wshape = 2, wscale = 1/gamma(1.5), cmu = 1, ctau = 1, sigmau = 1, xi = 0.1),
#'  xlim = c(0, 15), type = 'l', lwd = 2, 
#'   xlab = "x", ylab = "F(x)")
#' lines(xx, pgpd(xx, sigmau = 1, xi = 0.1), col = "red", lty = 2, lwd = 2)
#' lines(xx, pweibull(xx, shape = 2, scale = 1/gamma(1.5)), col = "blue", lty = 2, lwd = 2)
#' legend('bottomright', c('DWM', 'Weibull', 'GPD'),
#'       col = c("black", "blue", "red"), lty = c(1, 2, 2), lwd = 2)
#' }
#' 
NULL

#' @export
#' @aliases dwm ddwm pdwm qdwm rdwm
#' @rdname  dwm

# probability density function for dynamically weighted mixture model
ddwm = function(x, wshape = 1, wscale = 1, cmu = 1, ctau = 1,
  sigmau = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2),
  xi = 0, log = FALSE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.posparam(wshape, allowvec = TRUE)
  check.posparam(wscale, allowvec = TRUE)
  check.param(cmu, allowvec = TRUE) # not neccessarily positive
  check.posparam(ctau, allowvec = TRUE, allowzero = TRUE)
  check.posparam(sigmau, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.logic(log)

  n = check.inputn(c(length(x), length(wshape), length(wscale),
    length(cmu), length(ctau), length(sigmau), length(xi)), allowscalar = TRUE)
  oneparam = (check.inputn(c(length(wshape), length(wscale),
    length(cmu), length(ctau), length(sigmau), length(xi)), allowscalar = TRUE) == 1)
  
  if (any(is.infinite(x))) warning("infinite quantiles set to NA")

  x[is.infinite(x)] = NA # user will have to deal with infinite cases

  x = rep(x, length.out = n)
  wshape = rep(wshape, length.out = n)
  wscale = rep(wscale, length.out = n)
  cmu = rep(cmu, length.out = n)
  ctau = rep(ctau, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
    
  rx <- function(x, wshape, wscale, cmu, ctau, sigmau, xi) {
    (dgpd(x, 0, sigmau, xi) - dweibull(x, wshape, wscale))*atan((x - cmu)/ctau)
  }
  
  d = x # will pass through NA/NaN as entered

  whichnonmiss = which(!is.na(x))

  # As numerical integration is required for normalisation constant,
  # separate out case of scalar parameters in which this only needs to be calculated once
  if (oneparam) {
    r = try(integrate(rx, wshape = wshape[1], wscale = wscale[1],
      cmu = cmu[1], ctau = ctau[1], sigmau = sigmau[1], xi = xi[1],
      lower = 0, upper = Inf, subdivisions = 10000, rel.tol = 1e-10, stop.on.error = FALSE)$value)         
  
    if (inherits(r, "try-error")) {
      z = rep(NA, n)
    } else {              
      z = rep(1 + r/pi, length.out = n)
    }
  } else {
    z = rep(NA, n)
    for (i in 1:n) {
      r = try(integrate(r, wshape = wshape[i], wscale = wscale[i],
        cmu = cmu[i], ctau = ctau[i], sigmau = sigmau[i], xi = xi[i],
        lower = 0, upper = Inf, subdivisions = 10000, rel.tol = 1e-10, stop.on.error = FALSE)$value)         

      if (inherits(r, "try-error")) {
        z[i] = NA
      } else {              
        z[i] = 1 + r/pi
      }
    }
  }
  
  pweights = pcauchy(x[whichnonmiss], cmu[whichnonmiss], ctau[whichnonmiss])
  d[whichnonmiss] = ((1 - pweights) * dweibull(x[whichnonmiss], wshape[whichnonmiss], wscale[whichnonmiss]) +
    pweights * dgpd(x, 0, sigmau[whichnonmiss], xi[whichnonmiss]))/z[whichnonmiss]
  
  if (log) d = log(d)

  d
}

#' @export
#' @aliases dwm ddwm pdwm qdwm rdwm
#' @rdname  dwm

# cumulative distribution function for dynamically weighted mixture model
pdwm = function(q, wshape = 1, wscale = 1, cmu = 1, ctau = 1,
  sigmau = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2),
  xi = 0, lower.tail = TRUE) {
  
  # Check properties of inputs
  check.quant(q, allowna = TRUE, allowinf = TRUE)
  check.posparam(wshape, allowvec = TRUE)
  check.posparam(wscale, allowvec = TRUE)
  check.param(cmu, allowvec = TRUE) # not neccessarily positive
  check.posparam(ctau, allowvec = TRUE, allowzero = TRUE)
  check.posparam(sigmau, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.logic(lower.tail)

  n = check.inputn(c(length(q), length(wshape), length(wscale),
    length(cmu), length(ctau), length(sigmau), length(xi)), allowscalar = TRUE)
  oneparam = (check.inputn(c(length(wshape), length(wscale),
    length(cmu), length(ctau), length(sigmau), length(xi)), allowscalar = TRUE) == 1)

  if (any(is.infinite(q))) warning("infinite quantiles set to NA")

  q[is.infinite(q)] = NA # user will have to deal with infinite cases

  q = rep(q, length.out = n)
  wshape = rep(wshape, length.out = n)  
  wscale = rep(wscale, length.out = n)
  cmu = rep(cmu, length.out = n)
  ctau = rep(ctau, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
    
  rx <- function(x, wshape, wscale, cmu, ctau, sigmau, xi) {
    (dgpd(x, 0, sigmau, xi) - dweibull(x, wshape, wscale))*atan((x - cmu)/ctau)
  }

  rxw <- function(x, wshape, wscale, cmu, ctau) {
     (1 - pcauchy(x, cmu, ctau)) * dweibull(x, wshape, wscale)
  }

  rxg <- function(x, cmu, ctau, sigmau, xi){
     pcauchy(x, cmu, ctau) * dgpd(x, 0, sigmau, xi)
  }

  p = q # will pass through NA/NaN as entered
  
  whichnonmiss = which(!is.na(q))

  z1 = z2 = z = rep(NA, n)
  for (i in 1:n) {

    # As numerical integration is required for normalisation constant,
    # separate out case of scalar parameters in which this only needs to be calculated once
    if (oneparam & (i == 1)) {
      r = try(integrate(rx, wshape = wshape[i], wscale = wscale[i],
        cmu = cmu[i], ctau = ctau[i], sigmau = sigmau[i], xi = xi[i],
        lower = 0, upper = Inf, subdivisions = 10000, rel.tol = 1e-10, stop.on.error = FALSE)$value)         
  
      if (inherits(r, "try-error")) {
        z[i] = NA
      } else {              
        z[i] = 1 + r/pi
      }
    } else if (oneparam & (i > 1)) {
      z[i] = z[1]      
    } else {
      r = try(integrate(rx, wshape = wshape[i], wscale = wscale[i],
        cmu = cmu[i], ctau = ctau[i], sigmau = sigmau[i], xi = xi[i],
        lower = 0, upper = Inf, subdivisions = 10000, rel.tol = 1e-10, stop.on.error = FALSE)$value)         
  
      if (inherits(r, "try-error")) {
        z[i] = NA
      } else {              
        z[i] = 1 + r/pi
      }
    }
    
    r1 = try(integrate(rxw, wshape = wshape[i], wscale = wscale[i], cmu = cmu[i], ctau = ctau[i],
      lower = 0, upper = q[i], subdivisions = 10000, rel.tol = 1e-10, stop.on.error = FALSE)$value)         

    if (inherits(r1, "try-error")) {
      z1[i] = NA
    } else {              
      z1[i] = r1
    }

    r2 = try(integrate(rxg, cmu = cmu[i], ctau = ctau[i], sigmau = sigmau[i], xi = xi[i],
      lower = 0, upper = q[i], subdivisions = 10000, rel.tol = 1e-10, stop.on.error = FALSE)$value)         

    if (inherits(r2, "try-error")) {
      z2[i] = NA
    } else {              
      z2[i] = r2
    }
  }

  p[whichnonmiss] = (z1[whichnonmiss] + z2[whichnonmiss])/z[whichnonmiss]

  if (!lower.tail) p = 1 - p

  p
}

#' @export
#' @aliases dwm ddwm pdwm qdwm rdwm
#' @rdname  dwm

# inverse cumulative distribution function for dynamically weighted mixture model
qdwm = function(p, wshape = 1, wscale = 1, cmu = 1, ctau = 1,
  sigmau = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2),
  xi = 0, lower.tail = TRUE, qinit = NULL) {
  
  # Check properties of inputs
  check.prob(p, allowna = TRUE)
  check.posparam(wshape, allowvec = TRUE)
  check.posparam(wscale, allowvec = TRUE)
  check.param(cmu, allowvec = TRUE) # not neccessarily positive
  check.posparam(ctau, allowvec = TRUE, allowzero = TRUE)
  check.posparam(sigmau, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.logic(lower.tail)

  n = check.inputn(c(length(p), length(wshape), length(wscale),
    length(cmu), length(ctau), length(sigmau), length(xi)), allowscalar = TRUE)

  if (!lower.tail) p = 1 - p
  
  check.posparam(qinit, allowvec = TRUE, allownull = TRUE, allowzero = TRUE)
  if (is.null(qinit)) qinit = NA

  qinit = rep(qinit, length.out = n)
  
  p = rep(p, length.out = n)
  wshape = rep(wshape, length.out = n)
  wscale = rep(wscale, length.out = n)
  cmu = rep(cmu, length.out = n)
  ctau = rep(ctau, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  
  # No closed form solution for quantile function, need to solve numerically
  pdmmmin = function(q, cprob, wshape, wscale, cmu, ctau, sigmau, xi) {
     cdfmm = pdwm(q, wshape, wscale, cmu, ctau, sigmau, xi)
     if (is.na(cdfmm)) {
        qdiff = 1e6
     } else {
        qdiff = abs(cdfmm - cprob)
     }
    qdiff
  }
  
  findqdmm = function(cprob, wshape, wscale, cmu, ctau, sigmau, xi, qinit) {
    if (is.na(qinit)) {
      qwbl = qweibull(cprob, wshape, wscale)
      qgp = qgpd(cprob, 0, sigmau, xi)
      qinit = mean(c(qwbl, qgp))
    }
    
    gt = try(nlm(pdmmmin, qinit, cprob, wshape, wscale, cmu, ctau, sigmau, xi,
      gradtol = 1e-10, steptol = 1e-10)$estimate)

    if (inherits(gt, "try-error")) {
      gt = try(nlm(pdmmmin, qgpd(cprob, 0, sigmau, xi), cprob, wshape, wscale, cmu, ctau, sigmau, xi,
        gradtol = 1e-10, steptol = 1e-10)$estimate)
      
      if (inherits(gt, "try-error")) {
        gt = NA
      }
    }
    return(gt)
  }
 
  q = rep(NA, n)
  for (i in 1:n) {
    q[i] = findqdmm(p[i], wshape[i], wscale[i], cmu[i], ctau[i], sigmau[i], xi[i], qinit[i])
  }     

  q                     
}

#' @export
#' @aliases dwm ddwm pdwm qdwm rdwm
#' @rdname  dwm

# random number generation for dynamically weighted mixture model
rdwm = function(n = 1, wshape = 1, wscale = 1, cmu = 1, ctau = 1,
  sigmau = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2), xi = 0) {
  
  # Check properties of inputs
  check.n(n)
  check.posparam(wshape, allowvec = TRUE)
  check.posparam(wscale, allowvec = TRUE)
  check.param(cmu, allowvec = TRUE) # not neccessarily positive
  check.posparam(ctau, allowvec = TRUE, allowzero = TRUE)
  check.posparam(sigmau, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)

  check.inputn(c(n, length(wshape), length(wscale),
    length(cmu), length(ctau), length(sigmau), length(xi)), allowscalar = TRUE)

  if (any(xi == 1)) stop("shape cannot be 1")
  
  wshape = rep(wshape, length.out = n)
  wscale = rep(wscale, length.out = n)
  cmu = rep(cmu, length.out = n)
  ctau = rep(ctau, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  
  # Simulation scheme proposed by authors (see help for reference)
  r = rep(NA, n)

  i = 1
  while (i <= n) {
    u = runif(1)
    if (u < 0.5) {
      rw = rweibull(1, wshape[i], wscale[i])
      pw = pcauchy(rw, cmu[i], ctau[i])

      # accept or reject
      v = runif(1)
      if (v <= (1 - pw)) {
        r[i] = rw
        i = i + 1
      }
    } else {
      rg = rgpd(1, 0, sigmau[i], xi[i])
      pg = pcauchy(rg, cmu[i], ctau[i])

      # accept or reject
      v = runif(1)
      if (v <= pg){
        r[i] = rg
        i = i + 1
      }
    }
  } 
  
  r
}
