#' @name kfun
#' 
#' @param truncpoint   upper endpoint as standardised location \code{x/lambda} 
#' @inheritParams kernels
#' @inheritParams checking
#' 
#' @title Various subsidiary kernel function, conversion of bandwidths and evaluating certain
#' kernel integrals.
#'
#' @description Functions for checking the inputs to the kernel functions, evaluating 
#' integrals \eqn{\int u^l K*(u) du} for \eqn{l = 0, 1, 2} and conversion between the two bandwidth
#' definitions.
#' 
#' @details Various boundary correction methods require integral of (partial moments of)
#' kernel within the range of support, over the range \eqn{[-1, p]} where \eqn{p}
#' is the \code{truncpoint} determined by the standardised distance of location \eqn{x}
#' where KDE is being evaluated to the lower bound of zero, i.e. \code{truncpoint = x/lambda}.
#' The exception is the normal kernel which has unbounded support so the \eqn{[-5*\lambda, p]} where
#' \code{lambda} is the standard deviation bandwidth. There is a function for each partial moment
#' of degree (0, 1, 2):
#' \itemize{
#'  \item \code{ka0} - \eqn{\int_{-1}^{p} K*(z) dz}
#'  \item \code{ka1} - \eqn{\int_{-1}^{p} u K*(z) dz}
#'  \item \code{ka2} - \eqn{\int_{-1}^{p} u^2 K*(z) dz}
#' }
#' Notice that when evaluated at the upper endpoint on the support \eqn{p = 1}
#' (or \eqn{p = \infty} for normal) these are the zeroth, first and second moments. In the
#' normal distribution case the lower bound on the region of integration is \eqn{\infty} but
#' implemented here as \eqn{-5*\lambda}. 
#' These integrals are all specified in closed form, there is no need for numerical integration
#' (except normal which uses the \code{\link[stats:Normal]{pnorm}} function). 
#' 
#' See \code{\link[evmix:kernels]{kpu}} for list of kernels and discussion of bandwidth 
#' definitions (and their default values):
#' \enumerate{
#'  \item \code{bw} - in terms of number of standard deviations of the kernel, consistent
#'    with the defined values in the \code{\link[stats:density]{density}} function in
#'    the \code{R} base libraries
#'  \item \code{lambda} - in terms of half-width of kernel
#' }
#' The \code{\link[evmix:kfun]{klambda}} function converts the \code{bw} to the \code{lambda}
#' equivalent, and \code{\link[evmix:kfun]{kbw}} applies converse. These conversions are
#' kernel specific as they depend on the kernel standard deviations. If both \code{bw} and
#' \code{lambda} are provided then the latter is used by default. If neither are provided 
#' (\code{bw=NULL} and \code{lambda=NULL}) then default is \code{lambda=1}.
#' 
#' \code{\link[evmix:kfun]{check.kinputs}} checks all the kernel function inputs,
#' \code{\link[evmix:kfun]{check.klambda}} checks the pair of inputted bandwidths and
#' \code{\link[evmix:kfun]{check.kernel}} checks the kernel names.
#' 
#' @return \code{\link[evmix:kfun]{klambda}} and \code{\link[evmix:kfun]{kbw}} return the
#' \code{lambda} and \code{bw} bandwidths respectively.
#' 
#' The checking functions \code{\link[evmix:kfun]{check.kinputs}},
#' \code{\link[evmix:kfun]{check.klambda}} and \code{\link[evmix:kfun]{check.kernel}}
#' will stop on errors and return no value.
#' 
#' \code{\link[evmix:kfun]{ka0}}, \code{\link[evmix:kfun]{ka1}} and \code{\link[evmix:kfun]{ka2}}
#' return the partial moment integrals specified above. 
#'
#' @references
#' \url{http://en.wikipedia.org/wiki/Kernel_density_estimation}
#' 
#' \url{http://en.wikipedia.org/wiki/Kernel_(statistics)}
#' 
#' Wand and Jones (1995). Kernel Smoothing. Chapman & Hall.
#' 
#' @author Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}.
#'
#' @seealso \code{\link[evmix:kernels]{kernels}}, \code{\link[stats:density]{density}}, 
#' \code{\link[evmix:kden]{kden}} and \code{\link[evmix:bckden]{bckden}}.
#' 
#' @aliases klambda kbw check.kinputs check.kernel check.kbw ka0 ka1 ka2
#' @family kernels
#' 
#' @examples
#' xx = seq(-2, 2, 0.01)
#' plot(xx, kdgaussian(xx), type = "l", col = "black",ylim = c(0, 1.2))
#' lines(xx, kduniform(xx), col = "grey")
#' lines(xx, kdtriangular(xx), col = "blue")
#' lines(xx, kdepanechnikov(xx), col = "darkgreen")
#' lines(xx, kdbiweight(xx), col = "red")
#' lines(xx, kdtriweight(xx), col = "purple")
#' lines(xx, kdtricube(xx), col = "orange")
#' lines(xx, kdparzen(xx), col = "salmon")
#' lines(xx, kdcosine(xx), col = "cyan")
#' lines(xx, kdoptcosine(xx), col = "goldenrod")
#' legend("topright", c("Gaussian", "uniform", "triangular", "Epanechnikov",
#' "biweight", "triweight", "tricube", "Parzen", "cosine", "optcosine"), lty = 1,
#' col = c("black", "grey", "blue", "darkgreen", "red", "purple",
#'   "salmon", "orange", "cyan", "goldenrod"))
#' 
NULL

#' @export
#' @rdname kfun
check.kinputs <- function(x, lambda, bw, kerncentres, allownull = FALSE) {

  check.quant(x)
  check.quant(kerncentres)
  check.logic(allownull)

  if (all(c(length(x), length(kerncentres)) != 1))
    warnings("Check that both x and kerncentres are supposed to be vectors")
  
  check.kbw(lambda, bw, allownull = allownull)
}

#' @export
#' @rdname kfun
check.kernel <- function(kernel) {

  check.text(kernel)

  allkernels = c("gaussian", "normal", "uniform", "rectangular",
    "triangular", "epanechnikov", "biweight", "triweight", "tricube",
    "parzen", "cosine", "optcosine")
  
  if (!(kernel %in% allkernels))
    stop(paste(c("kernel must be one of", allkernels), collapse = " "))
}

#' @export
#' @rdname kfun
check.kbw <- function(lambda, bw, allownull = FALSE) {

  check.logic(allownull)

  if (!allownull & is.null(lambda) & is.null(bw)) stop("lambda and bw cannot both be NULL")
  
  if (!is.null(lambda)) {
    if (length(lambda) != 1) stop("lambda must be a scalar")
    if (mode(lambda) != "numeric") stop("lambda must be a scalar")
    if (!is.finite(lambda)) stop("lambda must be a finite scalar")
    if (lambda <= 0) stop("lambda must be a positive scalar")
  }
  
  if (!is.null(bw)) {
    if (length(bw) != 1) stop("bw must be a scalar")
    if (mode(bw) != "numeric") stop("bw must be a scalar")
    if (!is.finite(bw)) stop("bw must be a finite scalar")
    if (bw <= 0) stop("bw must be a positive scalar")
  }
}

#' @export
#' @rdname kfun
klambda <- function(bw = NULL, kernel = "gaussian", lambda = NULL) {
  # converts bw to equivalent lambda, and defaults to lambda if both given
  # lambda is limit of kernel for all except Gaussian
  # bw is in units of kernel standard deviations

  check.kernel(kernel)
  check.kbw(lambda, bw, allownull = TRUE) # allow both missing to default to lambda=1

  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  # variance of standard kernel
  kvar = switch(kernel,
    gaussian = 1,
    uniform = 1/3,
    triangular = 1/6,
    epanechnikov = 1/5,
    biweight = 1/7,
    triweight = 1/9,
    tricube = 35/243,
    parzen = 1/12,
    cosine = 1/3 - 2/pi/pi,
    optcosine = 1 - 8/pi/pi)
  
  if (!is.null(bw) & is.null(lambda)) {# bw provided only
    lambda = bw/sqrt(kvar)
  } else if (is.null(bw) & is.null(lambda)) {# neither provided
    lambda = 1 # default value
  } # if both then bw changed to lambda equivalent

  lambda
}

#' @export
#' @rdname kfun
kbw <- function(lambda = NULL, kernel = "gaussian", bw = NULL) {
  # converts lambda to equivalent bw, and defaults to lambda if both given
  # lambda is limit of kernel for all except Gaussian
  # bw is in units of kernel standard deviations

  check.kernel(kernel)
  check.kbw(lambda, bw, allownull = TRUE) # allow both missing to default to lambda=1

  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  # variance of standard kernel
  kvar = switch(kernel,
    gaussian = 1,
    uniform = 1/3,
    triangular = 1/6,
    epanechnikov = 1/5,
    biweight = 1/7,
    triweight = 1/9,
    tricube = 35/243,
    parzen = 1/12,
    cosine = 1/3 - 2/pi/pi,
    optcosine = 1 - 8/pi/pi)
  
  if (is.null(bw) & !is.null(lambda)) {
    bw = sqrt(kvar)*lambda
  } else if (!is.null(bw) & !is.null(lambda)) {
    bw = sqrt(kvar)*lambda # if both then bw changed to lambda equivalent
  } else if (is.null(bw) & is.null(lambda)) {
    bw = sqrt(kvar) # bw default (lambda = 1)
  }

  bw
}

#' @export
#' @rdname kfun
ka0 <- function(truncpoint, kernel = "gaussian") {

  check.kernel(kernel)
  
  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  switch(kernel,
    gaussian = kpgaussian(truncpoint),
    uniform = kpuniform(truncpoint),
    triangular = kptriangular(truncpoint),
    epanechnikov = kpepanechnikov(truncpoint),
    biweight = kpbiweight(truncpoint),
    triweight = kptriweight(truncpoint),
    tricube = kptricube(truncpoint),
    parzen = kpparzen(truncpoint),
    cosine = kpcosine(truncpoint),
    optcosine = kpoptcosine(truncpoint))
}

#' @export
#' @rdname kfun
ka1 <- function(truncpoint, kernel = "gaussian") {
  
  check.kernel(kernel)
  
  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  if (kernel != "gaussian") truncpoint = pmax(pmin(truncpoint, 1), -1)
  
  switch(kernel,
    gaussian = -dnorm(truncpoint),
    uniform = (truncpoint^2 - 1)/4,
    triangular = ifelse(truncpoint <= 0, (3*truncpoint^2 + 2*truncpoint^3 - 1)/6, (3*truncpoint^2 - 2*truncpoint^3 - 1)/6),
    epanechnikov = 3*(2*truncpoint^2 - truncpoint^4 - 1)/16,
    biweight = 5*(3*truncpoint^2 - 3*truncpoint^4 + truncpoint^6 - 1)/32,
    triweight = 35*(4*truncpoint^2 - 6*truncpoint^4 + 4*truncpoint^6 - truncpoint^8 - 1)/256,
    tricube = ifelse(truncpoint <= 0, 70*(truncpoint^2/2 + 3*truncpoint^5/5 + 3*truncpoint^8/8 + truncpoint^11/11 + 0.1 - 3/8 + 1/11)/81,
      70*(truncpoint^2/2 - 3*truncpoint^5/5 + 3*truncpoint^8/8 - truncpoint^11/11 + 0.1 - 3/8 + 1/11)/81),
    parzen = ifelse(truncpoint < -0.5, 32*truncpoint^5 + 120*truncpoint^4 + 160*truncpoint^3 + 80*truncpoint^2 - 8,
      ifelse((truncpoint >= -0.5) & (truncpoint < 0), -96*truncpoint^5 - 120*truncpoint^4 + 40*truncpoint^2 - 7,
      ifelse((truncpoint >= 0) & (truncpoint < 0.5), 96*truncpoint^5 - 120*truncpoint^4 + 40*truncpoint^2 - 7,
      -32*truncpoint^5 + 120*truncpoint^4 - 160*truncpoint^3 + 80*truncpoint^2 - 8)))/60,
    cosine = (truncpoint^2/4 - 0.25 + truncpoint*sin(pi*truncpoint)/2/pi + (cos(pi*truncpoint) + 1)/2/pi/pi),
    optcosine = (truncpoint*sin(pi*truncpoint/2) - 1)/2 + cos(pi*truncpoint/2)/pi)
}

#' @export
#' @rdname kfun
ka2 <- function(truncpoint, kernel = "gaussian") {
  
  check.kernel(kernel)
  
  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  if (kernel != "gaussian") truncpoint = pmax(pmin(truncpoint, 1), -1)

  switch(kernel,
    gaussian = (ka0(truncpoint) + truncpoint*ka1(truncpoint)),
    uniform = (truncpoint^3 + 1)/6,
    triangular = ifelse(truncpoint <= 0, (4*truncpoint^3 + 3*truncpoint^4 + 1)/12, (4*truncpoint^3 - 3*truncpoint^4 + 1)/12),
    epanechnikov = (5*truncpoint^3 - 3*truncpoint^5 + 2)/20,
    biweight = (5*truncpoint^3 - 6*truncpoint^5 + 15*truncpoint^7/7 + 8/7)/16,
    triweight = (35*truncpoint^3/3 - 21*truncpoint^5 + 15*truncpoint^7 - 35*truncpoint^9/9 + 16/9)/32,
    tricube = ifelse(truncpoint <= 0, 70*(truncpoint^3/3 + truncpoint^6/2 + truncpoint^9/3 + truncpoint^12/12 + 1/12)/81,
      70*(truncpoint^3/3 - truncpoint^6/2 + truncpoint^9/3 - truncpoint^12/12 + 1/12)/81),
    parzen = ifelse(truncpoint < -0.5, 2/45 + 8*truncpoint^3/9 + 2*truncpoint^4 + 8*truncpoint^5/5 + 4*truncpoint^6/9,
      ifelse((truncpoint >= -0.5) & (truncpoint < 0), 1/24 + 4*truncpoint^3/9 - 8*truncpoint^5/5 - 4*truncpoint^6/3,
      ifelse((truncpoint >= 0) & (truncpoint < 0.5), 1/24 + 4*truncpoint^3/9 - 8*truncpoint^5/5 + 4*truncpoint^6/3,
      7/180 + 8*truncpoint^3/9 - 2*truncpoint^4 + 8*truncpoint^5/5 - 4*truncpoint^6/9))),
    cosine = ((truncpoint^3 + 1)/6 + truncpoint^2*sin(pi*truncpoint)/2/pi + (truncpoint*cos(pi*truncpoint) - 1)/pi/pi - sin(pi*truncpoint)/pi/pi/pi),
    optcosine = (truncpoint^2*sin(pi*truncpoint/2) + 1)/2 + 2*truncpoint*cos(pi*truncpoint/2)/pi - 4*(sin(pi*truncpoint/2) + 1)/pi/pi)
}
