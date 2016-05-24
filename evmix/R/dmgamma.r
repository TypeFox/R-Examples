#' @name mgamma
#' 
#' @title Mixture of Gammas Distribution
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the mixture of gammas distribution. The parameters
#'   are the multiple gamma shapes \code{mgshape} scales \code{mgscale} and weights \code{mgweights}.
#'
#' @inheritParams gammagpd
#' @inheritParams gpd
#' @param mgshape     mgamma shape (positive) as list or vector
#' @param mgscale     mgamma scale (positive) as list or vector
#' @param mgweight    mgamma weights (positive) as list or vector (\code{NULL} for equi-weighted)
#' 
#' @details Distribution functions for weighted mixture of gammas.
#' 
#' Suppose there are \eqn{M>=1} gamma components in the mixture model. If you 
#' wish to have a single (scalar) value for each parameter within each of the
#' \eqn{M} components then these can be input as a vector of length \eqn{M}. If
#' you wish to input a vector of values for each parameter within each of the
#' \eqn{M} components, then they are input as a list with each entry the
#' parameter object for each component (which can either be a scalar or
#' vector as usual). No matter whether they are input as a vector or list there
#' must be \eqn{M} elements in \code{mgshape} and \code{mgscale}, one for each
#' gamma mixture component. Further, any vectors in the list of parameters must
#' of the same length of the \code{x, q, p} or equal to the sample size \code{n}, where
#' relevant.
#' 
#' If \code{mgweight=NULL} then equal weights for each component are assumed. Otherwise, 
#' \code{mgweight} must be a list of the same length as \code{mgshape} and 
#' \code{mgscale}, filled with positive values. In the latter case, the weights are rescaled
#' to sum to unity.
#' 
#' The gamma is defined on the non-negative reals. Though behaviour at zero depends on
#' the shape (\eqn{\alpha}):
#' \itemize{
#'  \item \eqn{f(0+)=\infty} for \eqn{0<\alpha<1};
#'  \item \eqn{f(0+)=1/\beta} for \eqn{\alpha=1} (exponential);
#'  \item \eqn{f(0+)=0} for \eqn{\alpha>1};
#' }
#' where \eqn{\beta} is the scale parameter.
#' 
#' @return \code{\link[evmix:mgamma]{dmgamma}} gives the density, 
#' \code{\link[evmix:mgamma]{pmgamma}} gives the cumulative distribution function,
#' \code{\link[evmix:mgamma]{qmgamma}} gives the quantile function and 
#' \code{\link[evmix:mgamma]{rmgamma}} gives a random sample.
#' 
#' @note All inputs are vectorised except \code{log} and \code{lower.tail}, and
#' the gamma mixture parameters can be vectorised within the list. The main
#' inputs (\code{x}, \code{p} or \code{q}) and parameters must be either a
#' scalar or a vector. If vectors are provided they must all be of the same
#' length, and the function will be evaluated for each element of vector. In
#' the case of \code{\link[evmix:mgamma]{rmgamma}} any input vector must be of
#' length \code{n}. The only exception is when the parameters are single scalar
#' values, input as vector of length \eqn{M}.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:mgamma]{rmgamma}} is 1.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x},
#' \code{p} and \code{q} are passed through as is and infinite values are set to
#' \code{NA}. None of these are not permitted for the parameters.
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' \url{http://www.math.canterbury.ac.nz/~c.scarrott/evmix}
#' 
#' \url{http://en.wikipedia.org/wiki/Gamma_distribution}
#' 
#' \url{http://en.wikipedia.org/wiki/Mixture_model}
#' 
#' McLachlan, G.J. and Peel, D. (2000). Finite Mixture Models. Wiley.
#' 
#' @author Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @section Acknowledgments: Thanks to Daniela Laas, University of St Gallen, Switzerland for reporting various bugs in these functions.
#'
#' @seealso \code{\link[evmix:gammagpd]{gammagpd}}, \code{\link[evmix:gpd]{gpd}}
#' and \code{\link[stats:GammaDist]{dgamma}}
#' 
#' @aliases mgamma dmgamma pmgamma qmgamma rmgamma
#' @family  mgamma fmgamma
#'          gammagpd gammagpdcon fgammagpd fgammagpdcon normgpd fnormgpd
#'          mgammagpd mgammagpdcon fmgammagpd fmgammagpdcon 
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(2, 1))
#' 
#' n = 1000
#' x = rmgamma(n, mgshape = c(1, 6), mgscale = c(1,2), mgweight = c(1, 2))
#' xx = seq(-1, 40, 0.01)
#' 
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 40))
#' lines(xx, dmgamma(xx, mgshape = c(1, 6), mgscale = c(1, 2), mgweight = c(1, 2)))
#' 
#' # By direct simulation
#' n1 = rbinom(1, n, 1/3) # sample size from population 1
#' x = c(rgamma(n1, shape = 1, scale = 1), rgamma(n - n1, shape = 6, scale = 2))
#' 
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 40))
#' lines(xx, dmgamma(xx, mgshape = c(1, 6), mgscale = c(1, 2), mgweight = c(1, 2)))
#' }
#'
NULL

#' @export
#' @aliases mgamma dmgamma pmgamma qmgamma rmgamma
#' @rdname  mgamma

# probability density function for mixture of gammas
dmgamma <- function(x, mgshape = 1, mgscale = 1,  mgweight = NULL, log = FALSE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)

  # user may try to input parameters as vectors if only scalar for each component
  if (!is.list(mgshape) & is.vector(mgshape)) {
    check.posparam(mgshape, allowvec = TRUE)
    mgshape = as.list(mgshape)    
  }
  if (!is.list(mgscale) & is.vector(mgscale)) {
    check.posparam(mgscale, allowvec = TRUE)
    mgscale = as.list(mgscale)    
  }
  
  if (!is.list(mgscale) | !is.list(mgshape))
    stop("gamma mixture parameters must be either a vector (if scalar parameters) or lists of same length (i.e. one object in list per component)")

  if (!is.null(mgweight)) {
    if (!is.list(mgweight) & is.vector(mgweight)) {
      check.posparam(mgweight, allowvec = TRUE)
      mgweight = as.list(mgweight)
    }
    if (!is.list(mgweight))
      stop("gamma mixture parameters must be either a vector (if scalar parameters) or lists of same length (i.e. one object in list per component)")
  }
  
  # How many components in mixture
  allM = c(length(mgshape), length(mgscale))
  if (!is.null(mgweight)) allM = c(allM, length(mgweight))
  
  M = unique(allM)
  if (length(M) != 1)
    stop("gamma mixture parameters must be either a vector (if scalar parameters) or lists of same length (i.e. one object in list per component)")
  
  if (is.null(mgweight)) mgweight = as.list(rep(1/M, M))

  linputs = length(x)
  for (i in 1:M) {
    check.posparam(mgshape[[i]], allowvec = TRUE)
    check.posparam(mgscale[[i]], allowvec = TRUE)
    check.posparam(mgweight[[i]], allowvec = TRUE)
    linputs = c(linputs, length(mgshape[[i]]), length(mgscale[[i]]), length(mgweight[[i]]))
  }
  
  check.logic(log)
  n = check.inputn(linputs, allowscalar = TRUE)

  if (any(is.infinite(x))) warning("infinite quantiles set to NA")

  x[is.infinite(x)] = NA # user will have to deal with infinite cases

  x = rep(x, length.out = n)
  for (i in 1:M) {
    mgshape[[i]] = rep(mgshape[[i]], length.out = n)
    mgscale[[i]] = rep(mgscale[[i]], length.out = n)
    mgweight[[i]] = rep(mgweight[[i]], length.out = n)
  }

  # Renormalise weights (deal with scalar or vector parameters separately)
  if (max(linputs[-1]) > 1) {
    mgweight = lapply(mgweight, FUN = function(x, y) x/y, y = rowSums(simplify2array(mgweight)))
  } else {
    mgweight = lapply(mgweight, FUN = function(x, y) x/y, y = sum(sapply(mgweight, FUN = function(x) x[1])))
  }

  d = x # will pass through NA/NaN as entered
    
  whichb = which(is.finite(x))
  nb = length(whichb)
  
  if (nb > 0) {
    d[whichb] = 0
    for (i in 1:M) {
      d[whichb] = d[whichb] + dgamma(x[whichb], shape = mgshape[[i]][whichb], 
        scale = mgscale[[i]][whichb]) * mgweight[[i]][whichb]
    }
  }
  
  if (log) d = log(d)

  d
}

#' @export
#' @aliases mgamma dmgamma pmgamma qmgamma rmgamma
#' @rdname  mgamma

# cumulative distribution function for mixture of gammas
pmgamma <- function(q, mgshape = 1, mgscale = 1,  mgweight = NULL, lower.tail = TRUE) {

  # Check properties of inputs
  check.quant(q, allowna = TRUE, allowinf = TRUE)

  # user may try to input parameters as vectors if only scalar for each component
  if (!is.list(mgshape) & is.vector(mgshape)) {
    check.posparam(mgshape, allowvec = TRUE)
    mgshape = as.list(mgshape)    
  }
  if (!is.list(mgscale) & is.vector(mgscale)) {
    check.posparam(mgscale, allowvec = TRUE)
    mgscale = as.list(mgscale)    
  }
  
  if (!is.list(mgscale) | !is.list(mgshape))
    stop("gamma mixture parameters must be either a vector (if scalar parameters) or lists of same length (i.e. one object in list per component)")

  if (!is.null(mgweight)) {
    if (!is.list(mgweight) & is.vector(mgweight)) {
      check.posparam(mgweight, allowvec = TRUE)
      mgweight = as.list(mgweight)    
    }
    if (!is.list(mgweight))
      stop("gamma mixture parameters must be either a vector (if scalar parameters) or lists of same length (i.e. one object in list per component)")
  }
  
  # How many components in mixture
  allM = c(length(mgshape), length(mgscale))
  if (!is.null(mgweight)) allM = c(allM, length(mgweight))
  
  M = unique(allM)
  if (length(M) != 1)
    stop("gamma mixture parameters must be either a vector (if scalar parameters) or lists of same length (i.e. one object in list per component)")
  
  if (is.null(mgweight)) mgweight = as.list(rep(1/M, M))

  linputs = length(q)
  for (i in 1:M) {
    check.posparam(mgshape[[i]], allowvec = TRUE)
    check.posparam(mgscale[[i]], allowvec = TRUE)
    check.posparam(mgweight[[i]], allowvec = TRUE)
    linputs = c(linputs, length(mgshape[[i]]), length(mgscale[[i]]), length(mgweight[[i]]))
  }
  
  check.logic(lower.tail)
  n = check.inputn(linputs, allowscalar = TRUE)

  if (any(is.infinite(q))) warning("infinite quantiles set to NA")

  q[is.infinite(q)] = NA # user will have to deal with infinite cases

  q = rep(q, length.out = n)
  for (i in 1:M) {
    mgshape[[i]] = rep(mgshape[[i]], length.out = n)
    mgscale[[i]] = rep(mgscale[[i]], length.out = n)
    mgweight[[i]] = rep(mgweight[[i]], length.out = n)
  }

  # Renormalise weights (deal with scalar or vector parameters separately)
  if (max(linputs[-1]) > 1) {
    mgweight = lapply(mgweight, FUN = function(x, y) x/y, y = rowSums(simplify2array(mgweight)))
  } else {
    mgweight = lapply(mgweight, FUN = function(x, y) x/y, y = sum(sapply(mgweight, FUN = function(x) x[1])))
  }
  
  p = q # will pass through NA/NaN as entered
  
  whichb = which(is.finite(q))
  nb = length(whichb)
  
  if (nb > 0) {
    p[whichb] = rep(0, nb)
    for (i in 1:M) {
      p[whichb] = p[whichb] + pgamma(q[whichb], shape = mgshape[[i]][whichb],
        scale = mgscale[[i]][whichb]) * mgweight[[i]][whichb]
    }
  }

  if (!lower.tail) p = 1 - p

  p
}

#' @export
#' @aliases mgamma dmgamma pmgamma qmgamma rmgamma
#' @rdname  mgamma

# inverse cumulative distribution function for mixture of gammas
qmgamma <- function(p, mgshape = 1, mgscale = 1,  mgweight = NULL, lower.tail = TRUE) {

  # Check properties of inputs
  check.prob(p, allowna = TRUE)

  # user may try to input parameters as vectors if only scalar for each component
  if (!is.list(mgshape) & is.vector(mgshape)) {
    check.posparam(mgshape, allowvec = TRUE)
    mgshape = as.list(mgshape)    
  }
  if (!is.list(mgscale) & is.vector(mgscale)) {
    check.posparam(mgscale, allowvec = TRUE)
    mgscale = as.list(mgscale)    
  }
  
  if (!is.list(mgscale) | !is.list(mgshape))
    stop("gamma mixture parameters must be either a vector (if scalar parameters) or lists of same length (i.e. one object in list per component)")

  if (!is.null(mgweight)) {
    if (!is.list(mgweight) & is.vector(mgweight)) {
      check.posparam(mgweight, allowvec = TRUE)
      mgweight = as.list(mgweight)    
    }
    if (!is.list(mgweight))
      stop("gamma mixture parameters must be either a vector (if scalar parameters) or lists of same length (i.e. one object in list per component)")
  }
  
  # How many components in mixture
  allM = c(length(mgshape), length(mgscale))
  if (!is.null(mgweight)) allM = c(allM, length(mgweight))
  
  M = unique(allM)
  if (length(M) != 1)
    stop("gamma mixture parameters must be either a vector (if scalar parameters) or lists of same length (i.e. one object in list per component)")
  
  if (is.null(mgweight)) mgweight = as.list(rep(1/M, M))

  linputs = length(p)
  for (i in 1:M) {
    check.posparam(mgshape[[i]], allowvec = TRUE)
    check.posparam(mgscale[[i]], allowvec = TRUE)
    check.posparam(mgweight[[i]], allowvec = TRUE)
    linputs = c(linputs, length(mgshape[[i]]), length(mgscale[[i]]), length(mgweight[[i]]))
  }
  
  check.logic(lower.tail)
  n = check.inputn(linputs, allowscalar = TRUE)

  if (!lower.tail) p = 1 - p

  p = rep(p, length.out = n)
  for (i in 1:M) {
    mgshape[[i]] = rep(mgshape[[i]], length.out = n)
    mgscale[[i]] = rep(mgscale[[i]], length.out = n)
    mgweight[[i]] = rep(mgweight[[i]], length.out = n)
  }

  # Renormalise weights (deal with scalar or vector parameters separately)
  if (max(linputs[-1]) > 1) {
    mgweight = lapply(mgweight, FUN = function(x, y) x/y, y = rowSums(simplify2array(mgweight)))
  } else {
    mgweight = lapply(mgweight, FUN = function(x, y) x/y, y = sum(sapply(mgweight, FUN = function(x) x[1])))
  }
  
  q = p # will pass through NA/NaN as entered
  
  whichb = which(is.finite(p))
  nb = length(whichb)

  # obtain quantile function below threshold by interpolation
  # - CDF is quick to calculate, interpolation algorithms also quite fast
  # - an alternative solution is to numerical solve to find quantile, but this is much slower
  # 
  # Treat single value of each parameter as special case as this much quicker,
  # parameter vectors have to be looped over which is painfully slow
  if (nb > 0) {
    maxp = max(ifelse(p[whichb] >= 1, NA, p[whichb]), na.rm = TRUE)

    # if only scalar parameters
    if (all(linputs[-1] == 1)) {
      maxq = max(qgamma(maxp, 
        shape = sapply(mgshape, FUN = function(x) x[1]),
        scale = sapply(mgscale, FUN = function(x) x[1])))
        
      qk = seq(0, maxq, length.out = 1000)

      pk = pmgamma(qk,
        sapply(mgshape, FUN = function(x) x[1]),
        sapply(mgscale, FUN = function(x) x[1]),
        sapply(mgweight, FUN = function(x) x[1])) 

      qfun = splinefun(x = pk, y = qk)

      q[whichb] = qfun(p[whichb])
    } else { # loop over each set of parameter values
      for (i in whichb) {
        maxq = max(qgamma(maxp,
          shape = sapply(mgshape, FUN = function(x, i) x[i], i = i),
          scale = sapply(mgscale, FUN = function(x, i) x[i], i = i)))
        
        qk = seq(0, maxq, length.out = 1000)
          
        pk = sapply(qk, FUN = pmgamma, 
          mgshape = sapply(mgshape, FUN = function(x, i) x[i], i = i),
          mgscale = sapply(mgscale, FUN = function(x, i) x[i], i = i),
          mgweight = sapply(mgweight, FUN = function(x, i) x[i], i = i)) 

        qfun = splinefun(x = pk, y = qk)

        q[i] = qfun(p[i])
      }
    }
  }

  q
}

#' @export
#' @aliases mgamma dmgamma pmgamma qmgamma rmgamma
#' @rdname  mgamma

# random number generation for mixture of gammas
rmgamma <- function(n = 1, mgshape = 1, mgscale = 1,  mgweight = NULL) {

  # Check properties of inputs
  check.n(n)

  # user may try to input parameters as vectors if only scalar for each component
  if (!is.list(mgshape) & is.vector(mgshape)) {
    check.posparam(mgshape, allowvec = TRUE)
    mgshape = as.list(mgshape)    
  }
  if (!is.list(mgscale) & is.vector(mgscale)) {
    check.posparam(mgscale, allowvec = TRUE)
    mgscale = as.list(mgscale)    
  }
  
  if (!is.list(mgscale) | !is.list(mgshape))
    stop("gamma mixture parameters must be either a vector (if scalar parameters) or lists of same length (i.e. one object in list per component)")

  if (!is.null(mgweight)) {
    if (!is.list(mgweight) & is.vector(mgweight)) {
      check.posparam(mgweight, allowvec = TRUE)
      mgweight = as.list(mgweight)    
    }
    if (!is.list(mgweight))
      stop("gamma mixture parameters must be either a vector (if scalar parameters) or lists of same length (i.e. one object in list per component)")
  }
  
  # How many components in mixture
  allM = c(length(mgshape), length(mgscale))
  if (!is.null(mgweight)) allM = c(allM, length(mgweight))
  
  M = unique(allM)
  if (length(M) != 1)
    stop("gamma mixture parameters must be either a vector (if scalar parameters) or lists of same length (i.e. one object in list per component)")
  
  if (is.null(mgweight)) mgweight = as.list(rep(1/M, M))

  linputs = n
  for (i in 1:M) {
    check.posparam(mgshape[[i]], allowvec = TRUE)
    check.posparam(mgscale[[i]], allowvec = TRUE)
    check.posparam(mgweight[[i]], allowvec = TRUE)
    linputs = c(linputs, length(mgshape[[i]]), length(mgscale[[i]]), length(mgweight[[i]]))
  }

  check.inputn(linputs, allowscalar = TRUE)

  for (i in 1:M) {
    mgshape[[i]] = rep(mgshape[[i]], length.out = n)
    mgscale[[i]] = rep(mgscale[[i]], length.out = n)
    mgweight[[i]] = rep(mgweight[[i]], length.out = n)
  }

  # Renormalise weights (deal with scalar or vector parameters separately)
  if (max(linputs[-1]) > 1) {
    mgweight = lapply(mgweight, FUN = function(x, y) x/y, y = rowSums(simplify2array(mgweight)))
  } else {
    mgweight = lapply(mgweight, FUN = function(x, y) x/y, y = sum(sapply(mgweight, FUN = function(x) x[1])))
  }
  
  # Treat scalar and vector parameters seperately
  # if only scalar parameters
  if (all(linputs[-1] == 1)) {
    nM = rmultinom(1, n, prob = sapply(mgweight, FUN = function(x) x[1]))
    
    shapeM = rep(sapply(mgshape, FUN = function(x) x[1]), times = nM)
    scaleM = rep(sapply(mgscale, FUN = function(x) x[1]), times = nM)

    sample(rgamma(n, shape = shapeM, scale = scaleM))
  } else {
    
    allweights = simplify2array(mgweight)
    whichM = apply(allweights, MARGIN = 1, 
      FUN = function(weights, M) sample(1:M, 1, prob = weights), M = M)
    
    whichparam = which(col(allweights) == whichM)
    
    shapeM = simplify2array(mgshape)[whichparam]
    scaleM = simplify2array(mgscale)[whichparam]
    
    sample(rgamma(n, shape = shapeM, scale = scaleM))
  }
}
