#' @name mgammagpdcon
#' 
#' @title Mixture of Gammas Bulk and GPD Tail Extreme Value Mixture Model with Single Continuity Constraint
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the extreme value mixture model with mixture of gammas for bulk
#'   distribution upto the threshold and conditional GPD for upper tail with continuity at threshold. The parameters
#'   are the multiple gamma shapes \code{mgshape}, scales \code{mgscale} and \code{mgweights}, threshold \code{u}
#'   GPD shape \code{xi} and tail fraction \code{phiu}.
#'
#' @inheritParams mgammagpd
#' 
#' @details Extreme value mixture model combining mixture of gammas for the bulk
#' below the threshold and GPD for upper tail with continuity at threshold. 
#' 
#' The user can pre-specify \code{phiu} permitting a parameterised value for the tail
#' fraction \eqn{\phi_u}. Alternatively, when \code{phiu=TRUE} the tail fraction is
#' estimated as the tail fraction from the mixture of gammas bulk model.
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
#' The cumulative distribution function with tail fraction \eqn{\phi_u} defined by the
#' upper tail fraction of the mixture of gammas bulk model (\code{phiu=TRUE}), upto the 
#' threshold \eqn{0 < x \le u}, given by:
#' \deqn{F(x) = H(x)}
#' and above the threshold \eqn{x > u}:
#' \deqn{F(x) = H(u) + [1 - H(u)] G(x)}
#' where \eqn{H(x)} and \eqn{G(X)} are the mixture of gammas and conditional GPD
#' cumulative distribution functions.
#' 
#' The cumulative distribution function for pre-specified \eqn{\phi_u}, upto the
#' threshold \eqn{0 < x \le u}, is given by:
#' \deqn{F(x) = (1 - \phi_u) H(x)/H(u)}
#' and above the threshold \eqn{x > u}:
#' \deqn{F(x) = \phi_u + [1 - \phi_u] G(x)}
#' Notice that these definitions are equivalent when \eqn{\phi_u = 1 - H(u)}.
#' 
#' The continuity constraint means that \eqn{(1 - \phi_u) h(u)/H(u) = \phi_u g(u)}
#' where \eqn{h(x)} and \eqn{g(x)} are the mixture of gammas and conditional GPD
#' density functions respectively. The resulting GPD scale parameter is then:
#' \deqn{\sigma_u = \phi_u H(u) / [1 - \phi_u] h(u)}.
#' In the special case of where the tail fraction is defined by the bulk model this reduces to
#' \deqn{\sigma_u = [1 - H(u)] / h(u)}.
#' 
#' The gamma is defined on the non-negative reals, so the threshold must be positive. 
#' Though behaviour at zero depends on the shape (\eqn{\alpha}):
#' \itemize{
#'  \item \eqn{f(0+)=\infty} for \eqn{0<\alpha<1};
#'  \item \eqn{f(0+)=1/\beta} for \eqn{\alpha=1} (exponential);
#'  \item \eqn{f(0+)=0} for \eqn{\alpha>1};
#' }
#' where \eqn{\beta} is the scale parameter.
#' 
#' See \code{\link[evmix:gammagpd]{gammagpd}} for details of simpler parametric mixture model
#' with single gamma for bulk component and GPD for upper tail.
#' 
#' @return \code{\link[evmix:mgammagpdcon]{dmgammagpdcon}} gives the density, 
#' \code{\link[evmix:mgammagpdcon]{pmgammagpdcon}} gives the cumulative distribution function,
#' \code{\link[evmix:mgammagpdcon]{qmgammagpdcon}} gives the quantile function and 
#' \code{\link[evmix:mgammagpdcon]{rmgammagpdcon}} gives a random sample.
#' 
#' @note All inputs are vectorised except \code{log} and \code{lower.tail}, and the gamma mixture
#' parameters can be vectorised within the list. The main inputs (\code{x}, \code{p} or \code{q})
#' and parameters must be either a scalar or a vector. If vectors are provided they must all be
#' of the same length, and the function will be evaluated for each element of vector. In the case of 
#' \code{\link[evmix:mgammagpdcon]{rmgammagpdcon}} any input vector must be of length \code{n}.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:mgammagpdcon]{rmgammagpdcon}} is 1.
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
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' \url{http://en.wikipedia.org/wiki/Mixture_model}
#' 
#' McLachlan, G.J. and Peel, D. (2000). Finite Mixture Models. Wiley.
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' do Nascimento, F.F., Gamerman, D. and Lopes, H.F. (2011). A semiparametric
#' Bayesian approach to extreme value estimation. Statistical Computing, 22(2), 661-675.
#' 
#' @author Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @section Acknowledgments: Thanks to Daniela Laas, University of St Gallen, Switzerland for reporting various bugs in these functions.
#' 
#' @seealso \code{\link[evmix:gpd]{gpd}} and \code{\link[stats:GammaDist]{dgamma}}
#' 
#' @aliases mgammagpdcon dmgammagpdcon pmgammagpdcon qmgammagpdcon rmgammagpdcon
#' @family  mgamma fmgamma
#'          gammagpd gammagpdcon fgammagpd fgammagpdcon normgpd fnormgpd
#'          mgammagpd mgammagpdcon fmgammagpd fmgammagpdcon 
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(1, 1))
#' 
#' x = rmgammagpdcon(1000, mgshape = c(1, 6), mgscale = c(1, 2), mgweight = c(1, 2), u = 15, xi = 0)
#' xx = seq(-1, 40, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 40))
#' lines(xx, dmgammagpdcon(xx, mgshape = c(1, 6), mgscale = c(1, 2), mgweight = c(1, 2),
#'  u = 15, xi = 0))
#' abline(v = 15)
#' }
#'
NULL

#' @export
#' @aliases mgammagpdcon dmgammagpdcon pmgammagpdcon qmgammagpdcon rmgammagpdcon
#' @rdname  mgammagpdcon

# probability density function for mixture of gammas bulk with GPD for upper tail
dmgammagpdcon <- function(x, mgshape = 1, mgscale = 1,  mgweight = NULL,
  u = qgamma(0.9, mgshape[[1]], 1/mgscale[[1]]), xi = 0, phiu = TRUE, log = FALSE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.posparam(u, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.phiu(phiu, allowvec = TRUE)

  # user may try to input parameters as vectors if only scalar for each component
  if (!is.list(mgshape) & is.vector(mgshape)) {
    check.posparam(mgshape, allowvec = TRUE)
    mgshape = as.list(mgshape)    
  }
  if (!is.list(mgscale) & is.vector(mgscale)) {
    check.posparam(mgscale, allowvec = TRUE)
    mgscale = as.list(mgscale)    
  }
  
  if (!is.list(mgshape) | !is.list(mgshape))
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

  linputs = c(length(x), length(u), length(xi), length(phiu))
  for (i in 1:M) {
    check.posparam(mgshape[[i]], allowvec = TRUE)
    check.posparam(mgscale[[i]], allowvec = TRUE)
    check.posparam(mgweight[[i]], allowvec = TRUE)
    linputs = c(linputs, length(mgshape[[i]]), length(mgscale[[i]]), length(mgweight[[i]]))
  }
  
  check.logic(log)

  n = check.inputn(linputs, allowscalar = TRUE)

  # vectors of parameters or scalars? Latter simplifies calculations
  scalar.params = all(linputs[-1] == 1)
  
  if (any(is.infinite(x))) warning("infinite quantiles set to NA")

  x[is.infinite(x)] = NA # user will have to deal with infinite cases

  x = rep(x, length.out = n)
  for (i in 1:M) {
    mgshape[[i]] = rep(mgshape[[i]], length.out = n)
    mgscale[[i]] = rep(mgscale[[i]], length.out = n)
    mgweight[[i]] = rep(mgweight[[i]], length.out = n)
  }
  u = rep(u, length.out = n)
  xi = rep(xi, length.out = n)

  # Renormalise weights
  if (scalar.params) {
    mgweight = lapply(mgweight, FUN = function(x, y) x/y,
      y = sum(sapply(mgweight, FUN = function(x) x[1])))
    
    # Calculate CDF of mixture of gammas upto threshold (to get phiu)
    pu = pmgamma(u[1], sapply(mgshape, FUN = function(x) x[1]),
      sapply(mgscale, FUN = function(x) x[1]), sapply(mgweight, FUN = function(x) x[1]))
    du = dmgamma(u[1], sapply(mgshape, FUN = function(x) x[1]),
      sapply(mgscale, FUN = function(x) x[1]), sapply(mgweight, FUN = function(x) x[1]))
    
  } else {
    mgweight = lapply(mgweight, FUN = function(x, y) x/y, y = rowSums(simplify2array(mgweight)))

    # Calculate CDF of mixture of gammas upto threshold (to get phiu)
    pu = pmgamma(u, mgshape, mgscale, mgweight)
    du = dmgamma(u, mgshape, mgscale, mgweight)
  }

  if (is.logical(phiu)) {
    phiu = 1 - pu
  } else {
    phiu = rep(phiu, length.out = n)
  }
  phib = (1 - phiu) / pu

  sigmau = phiu / (phib * du)
  
  check.posparam(sigmau, allowvec = TRUE)

  if (scalar.params) {
    dmgammagpd(x, sapply(mgshape, FUN = function(x) x[1]), sapply(mgscale, FUN = function(x) x[1]),
      sapply(mgweight, FUN = function(x) x[1]), u[1], sigmau[1], xi[1], phiu[1], log)
  } else {
    dmgammagpd(x, mgshape, mgscale, mgweight, u, sigmau, xi, phiu, log)    
  }
}

#' @export
#' @aliases mgammagpdcon dmgammagpdcon pmgammagpdcon qmgammagpdcon rmgammagpdcon
#' @rdname  mgammagpdcon

# cumulative distribution function for mixture of gammas bulk with GPD for upper tail
pmgammagpdcon <- function(q, mgshape = 1, mgscale = 1,  mgweight = NULL,
  u = qgamma(0.9, mgshape[[1]], 1/mgscale[[1]]), xi = 0, phiu = TRUE, lower.tail = TRUE) {

  # Check properties of inputs
  check.quant(q, allowna = TRUE, allowinf = TRUE)
  check.posparam(u, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.phiu(phiu, allowvec = TRUE)

  # user may try to input parameters as vectors if only scalar for each component
  if (!is.list(mgshape) & is.vector(mgshape)) {
    check.posparam(mgshape, allowvec = TRUE)
    mgshape = as.list(mgshape)    
  }
  if (!is.list(mgscale) & is.vector(mgscale)) {
    check.posparam(mgscale, allowvec = TRUE)
    mgscale = as.list(mgscale)    
  }
  
  if (!is.list(mgshape) | !is.list(mgshape))
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

  linputs = c(length(q), length(u), length(xi), length(phiu))
  for (i in 1:M) {
    check.posparam(mgshape[[i]], allowvec = TRUE)
    check.posparam(mgscale[[i]], allowvec = TRUE)
    check.posparam(mgweight[[i]], allowvec = TRUE)
    linputs = c(linputs, length(mgshape[[i]]), length(mgscale[[i]]), length(mgweight[[i]]))
  }
  
  check.logic(lower.tail)

  n = check.inputn(linputs, allowscalar = TRUE)

  # vectors of parameters or scalars? Latter simplifies calculations
  scalar.params = all(linputs[-1] == 1)
  
  if (any(is.infinite(q))) warning("infinite quantiles set to NA")

  q[is.infinite(q)] = NA # user will have to deal with infinite cases

  q = rep(q, length.out = n)
  for (i in 1:M) {
    mgshape[[i]] = rep(mgshape[[i]], length.out = n)
    mgscale[[i]] = rep(mgscale[[i]], length.out = n)
    mgweight[[i]] = rep(mgweight[[i]], length.out = n)
  }
  u = rep(u, length.out = n)
  xi = rep(xi, length.out = n)

  # Renormalise weights
  if (scalar.params) {
    mgweight = lapply(mgweight, FUN = function(x, y) x/y,
      y = sum(sapply(mgweight, FUN = function(x) x[1])))
    
    # Calculate CDF of mixture of gammas upto threshold (to get phiu)
    pu = pmgamma(u[1], sapply(mgshape, FUN = function(x) x[1]),
      sapply(mgscale, FUN = function(x) x[1]), sapply(mgweight, FUN = function(x) x[1]))
    du = dmgamma(u[1], sapply(mgshape, FUN = function(x) x[1]),
      sapply(mgscale, FUN = function(x) x[1]), sapply(mgweight, FUN = function(x) x[1]))
    
  } else {
    mgweight = lapply(mgweight, FUN = function(x, y) x/y, y = rowSums(simplify2array(mgweight)))

    # Calculate CDF of mixture of gammas upto threshold (to get phiu)
    pu = pmgamma(u, mgshape, mgscale, mgweight)
    du = dmgamma(u, mgshape, mgscale, mgweight)
  }
  
  if (is.logical(phiu)) {
    phiu = 1 - pu
  } else {
    phiu = rep(phiu, length.out = n)
  }
  phib = (1 - phiu) / pu

  sigmau = phiu / (phib * du)
  
  check.posparam(sigmau, allowvec = TRUE)

  if (scalar.params) {
    pmgammagpd(q, sapply(mgshape, FUN = function(x) x[1]), sapply(mgscale, FUN = function(x) x[1]),
      sapply(mgweight, FUN = function(x) x[1]), u[1], sigmau[1], xi[1], phiu[1], lower.tail)
  } else {
    pmgammagpd(q, mgshape, mgscale, mgweight, u, sigmau, xi, phiu, lower.tail)    
  }
}

#' @export
#' @aliases mgammagpdcon dmgammagpdcon pmgammagpdcon qmgammagpdcon rmgammagpdcon
#' @rdname  mgammagpdcon

# inverse cumulative distribution function for mixture of gammas bulk with GPD for upper tail
qmgammagpdcon <- function(p, mgshape = 1, mgscale = 1, mgweight = NULL,
  u = qgamma(0.9, mgshape[[1]], 1/mgscale[[1]]), xi = 0, phiu = TRUE, lower.tail = TRUE) {

  # Check properties of inputs
  check.prob(p, allowna = TRUE)
  check.posparam(u, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.phiu(phiu, allowvec = TRUE)

  # user may try to input parameters as vectors if only scalar for each component
  if (!is.list(mgshape) & is.vector(mgshape)) {
    check.posparam(mgshape, allowvec = TRUE)
    mgshape = as.list(mgshape)    
  }
  if (!is.list(mgscale) & is.vector(mgscale)) {
    check.posparam(mgscale, allowvec = TRUE)
    mgscale = as.list(mgscale)    
  }
  
  if (!is.list(mgshape) | !is.list(mgshape))
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

  linputs = c(length(p), length(u), length(xi), length(phiu))
  for (i in 1:M) {
    check.posparam(mgshape[[i]], allowvec = TRUE)
    check.posparam(mgscale[[i]], allowvec = TRUE)
    check.posparam(mgweight[[i]], allowvec = TRUE)
    linputs = c(linputs, length(mgshape[[i]]), length(mgscale[[i]]), length(mgweight[[i]]))
  }
  
  check.logic(lower.tail)
  
  n = check.inputn(linputs, allowscalar = TRUE)

  if (!lower.tail) p = 1 - p

  # vectors of parameters or scalars? Latter simplifies calculations
  scalar.params = all(linputs[-1] == 1)
  
  p = rep(p, length.out = n)
  for (i in 1:M) {
    mgshape[[i]] = rep(mgshape[[i]], length.out = n)
    mgscale[[i]] = rep(mgscale[[i]], length.out = n)
    mgweight[[i]] = rep(mgweight[[i]], length.out = n)
  }
  u = rep(u, length.out = n)
  xi = rep(xi, length.out = n)

  # Renormalise weights
  if (scalar.params) {
    mgweight = lapply(mgweight, FUN = function(x, y) x/y,
      y = sum(sapply(mgweight, FUN = function(x) x[1])))
    
    # Calculate CDF of mixture of gammas upto threshold (to get phiu)
    pu = pmgamma(u[1], sapply(mgshape, FUN = function(x) x[1]),
      sapply(mgscale, FUN = function(x) x[1]), sapply(mgweight, FUN = function(x) x[1]))
    du = dmgamma(u[1], sapply(mgshape, FUN = function(x) x[1]),
      sapply(mgscale, FUN = function(x) x[1]), sapply(mgweight, FUN = function(x) x[1]))
    
  } else {
    mgweight = lapply(mgweight, FUN = function(x, y) x/y, y = rowSums(simplify2array(mgweight)))

    # Calculate CDF of mixture of gammas upto threshold (to get phiu)
    pu = pmgamma(u, mgshape, mgscale, mgweight)
    du = dmgamma(u, mgshape, mgscale, mgweight)
  }
  
  if (is.logical(phiu)) {
    phiu = 1 - pu
  } else {
    phiu = rep(phiu, length.out = n)
  }
  phib = (1 - phiu) / pu
  
  sigmau = phiu / (phib * du)
  
  check.posparam(sigmau, allowvec = TRUE)

  # for efficiency only want to pass scalar parameters, if they are input as scalar
  if (scalar.params) {
    qmgammagpd(p, sapply(mgshape, FUN = function(x) x[1]), sapply(mgscale, FUN = function(x) x[1]),
      sapply(mgweight, FUN = function(x) x[1]), u[1], sigmau[1], xi[1], phiu[1])
  } else {
    qmgammagpd(p, mgshape, mgscale, mgweight, u, sigmau, xi, phiu)
  }
}

#' @export
#' @aliases mgammagpdcon dmgammagpdcon pmgammagpdcon qmgammagpdcon rmgammagpdcon
#' @rdname  mgammagpdcon

# random number generation for mixture of gammas bulk with GPD for upper tail
rmgammagpdcon <- function(n = 1, mgshape = 1, mgscale = 1, mgweight = NULL,
  u = qgamma(0.9, mgshape[[1]], 1/mgscale[[1]]), xi = 0, phiu = TRUE) {

  # Check properties of inputs
  check.n(n)
  check.posparam(u, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.phiu(phiu, allowvec = TRUE)

  # user may try to input parameters as vectors if only scalar for each component
  if (!is.list(mgshape) & is.vector(mgshape)) {
    check.posparam(mgshape, allowvec = TRUE)
    mgshape = as.list(mgshape)    
  }
  if (!is.list(mgscale) & is.vector(mgscale)) {
    check.posparam(mgscale, allowvec = TRUE)
    mgscale = as.list(mgscale)    
  }
  
  if (!is.list(mgshape) | !is.list(mgshape))
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
  
  if (any(xi == 1)) stop("shape cannot be 1")

  qmgammagpdcon(runif(n), mgshape, mgscale, mgweight, u, xi, phiu)
}
