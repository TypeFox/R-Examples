#' @name kernels
#' 
#' @param x            location to evaluate KDE (single scalar or vector)
#' @param kerncentres  kernel centres (typically sample data vector or scalar)
#' @param lambda       bandwidth for kernel (as half-width of kernel) or \code{NULL}
#' @param bw           bandwidth for kernel (as standard deviations of kernel) or \code{NULL}
#' @param z            standardised location put into kernel \code{z = (x-kerncentres)/lambda}
#' @param kernel       kernel name (\code{default = "gaussian"})
#' 
#' @title Kernel functions
#'
#' @description Functions for commonly used kernels for kernel density estimation. The
#' density and cumulative distribution functions are provided.
#' 
#' @details Functions for the commonly used kernels for kernel density estimation. The
#' density and cumulative distribution functions are provided. Each function can accept the
#' bandwidth specified as either:
#' \enumerate{
#'  \item \code{bw} - in terms of number of standard deviations of the kernel, consistent
#'    with the defined values in the \code{\link[stats:density]{density}} function in
#'    the \code{R} base libraries
#'  \item \code{lambda} - in terms of half-width of kernel
#' }
#' If both bandwidths are given as \code{NULL} then the default bandwidth is \code{lambda=1}. If
#' either one is specified then this will be used. If both are specified then \code{lambda}
#' will be used.
#' 
#' All the kernels have bounded support \eqn{[-\lambda, \lambda]}, except the normal
#' (\code{"gaussian"}) which is unbounded. In the latter, both bandwidths are the same
#' \code{bw=lambda} and equal to the standard deviation.
#' 
#' Typically,a single location \code{x} at which to evaluate kernel is given along with
#' vector of kernel centres. As such, they are designed to be used with 
#' \code{\link[base:lapply]{sapply}} to loop over vector of locations at which to evaluate KDE. 
#' Alternatively, a vector of locations \code{x} can be given with a single scalar kernel centre
#' \code{kerncentres}, which is commonly used when locations are pre-standardised by
#' \code{(x-kerncentres)/lambda} and \code{kerncentre=0}. A warnings is given if both the
#' evaluation locations and kernel centres are vectors as this is not often needed so is
#' likely to be a user error.
#' 
#' If no kernel centres are provided then by default it is set to zero (i.e. x is at middle of kernel).
#' 
#' The following kernels are implemented, with relevant ones having definitions
#' consistent with those of the \code{\link[stats:density]{density}} function,
#' except where specified:
#' \itemize{
#'  \item \code{gaussian} or \code{normal}
#'  \item \code{uniform} or \code{rectangular} - same as \code{"rectangular"} in 
#'    \code{\link[stats:density]{density}} function
#'  \item \code{triangular}
#'  \item \code{epanechnikov}
#'  \item \code{biweight}
#'  \item \code{triweight}
#'  \item \code{tricube}
#'  \item \code{parzen}
#'  \item \code{cosine}
#'  \item \code{optcosine}
#' }
#' The kernel densities are all normalised to unity. See Wikipedia reference below
#' for their definitions.
#' 
#' Each kernel's functions can be called individually, or the global functions
#' \code{\link[evmix:kernels]{kdz}} and \code{\link[evmix:kernels]{kpz}} for the density and
#' cumulative distribution function can apply any particular kernel which is specified by the
#' \code{kernel} input. These global functions take the standardised locations
#' \code{z = (x - kerncentres)/lambda}.
#' 
#' @return code{\link[evmix:kernels]{kd*}}  and \code{\link[evmix:kernels]{kp*}} give the
#'  density and cumulative distribution functions for each kernel respectively, where
#'  \code{*} is the kernel name. \code{\link[evmix:kernels]{kdz}} and
#'  \code{\link[evmix:kernels]{kpz}} are the equivalent global functions for all of the 
#'  kernels.
#'  
#' @references
#' \url{http://en.wikipedia.org/wiki/Kernel_density_estimation}
#' 
#' \url{http://en.wikipedia.org/wiki/Kernel_(statistics)}
#' 
#' Wand, M. and Jones, M.C. (1995). Kernel Smoothing. Chapman && Hall.
#' 
#' @author Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}.
#'
#' @seealso \code{\link[stats:density]{density}}, \code{\link[evmix:kden]{kden}}
#' and \code{\link[evmix:bckden]{bckden}}.
#' 
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
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
#' col = c("black", "grey", "blue", "darkgreen", "red", "purple", "orange",
#'   "salmon", "cyan", "goldenrod"))
#' 
NULL

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kdgaussian <- function(x = 0, lambda = NULL, bw = NULL, kerncentres = 0) {
  
  check.kinputs(x, lambda, bw, kerncentres, allownull = TRUE)

  lambda = klambda(bw, "gaussian", lambda)

  z = (x - kerncentres)/lambda
  dnorm(z)/lambda
}

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kduniform <- function(x = 0, lambda = NULL, bw = NULL, kerncentres = 0) {
  
  check.kinputs(x, lambda, bw, kerncentres, allownull = TRUE)
  
  lambda = klambda(bw, "uniform", lambda)
  
  z = (x - kerncentres)/lambda
  (abs(z) <= 1)/lambda/2
}

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kdtriangular <- function(x = 0, lambda = NULL, bw = NULL, kerncentres = 0) {
  
  check.kinputs(x, lambda, bw, kerncentres, allownull = TRUE)
  
  lambda = klambda(bw, "triangular", lambda)
  
  z = (x - kerncentres)/lambda
  pmax(1 - abs(z), 0)/lambda
}

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kdepanechnikov <- function(x = 0, lambda = NULL, bw = NULL, kerncentres = 0) {
  
  check.kinputs(x, lambda, bw, kerncentres, allownull = TRUE)
  
  lambda = klambda(bw, "epanechnikov", lambda)
  
  z = (x - kerncentres)/lambda
  3*pmax(1 - z^2, 0)/4/lambda
}

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kdbiweight <- function(x = 0, lambda = NULL, bw = NULL, kerncentres = 0) {
  
  check.kinputs(x, lambda, bw, kerncentres, allownull = TRUE)
  
  lambda = klambda(bw, "biweight", lambda)
  
  z = (x - kerncentres)/lambda
  15*pmax(1 - z^2, 0)^2/16/lambda
}

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kdtriweight <- function(x = 0, lambda = NULL, bw = NULL, kerncentres = 0) {
  
  check.kinputs(x, lambda, bw, kerncentres, allownull = TRUE)
  
  lambda = klambda(bw, "triweight", lambda)
  
  z = (x - kerncentres)/lambda
  35*pmax(1 - z^2, 0)^3/32/lambda
}

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kdtricube <- function(x = 0, lambda = NULL, bw = NULL, kerncentres = 0) {
  
  check.kinputs(x, lambda, bw, kerncentres, allownull = TRUE)
  
  lambda = klambda(bw, "tricube", lambda)
  
  z = (x - kerncentres)/lambda
  70*pmax(1 - abs(z)^3, 0)^3/81/lambda
}

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kdparzen <- function(x = 0, lambda = NULL, bw = NULL, kerncentres = 0) {
  
  check.kinputs(x, lambda, bw, kerncentres, allownull = TRUE)
  
  lambda = klambda(bw, "parzen", lambda)
  
  z = (x - kerncentres)/lambda
  z = pmin(abs(z), 1)
  ((z > 0.5)*8*(1 - z)^3/3 + (z <= 0.5)*(4 - 24*z^2 + 24*z^3)/3)/lambda
}

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kdcosine <- function(x = 0, lambda = NULL, bw = NULL, kerncentres = 0) {
  
  check.kinputs(x, lambda, bw, kerncentres, allownull = TRUE)
  
  lambda = klambda(bw, "cosine", lambda)
  
  z = (x - kerncentres)/lambda
  (abs(z) <= 1)*(1 + cos(pi*z))/2/lambda
}

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kdoptcosine <- function(x = 0, lambda = NULL, bw = NULL, kerncentres = 0) {
  
  check.kinputs(x, lambda, bw, kerncentres, allownull = TRUE)
  
  lambda = klambda(bw, "optcosine", lambda)

  z = (x - kerncentres)/lambda
  pi*(abs(z) <= 1)*cos(pi*z/2)/4/lambda
}

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kpgaussian <- function(x = 0, lambda = NULL, bw = NULL, kerncentres = 0) {
  
  check.kinputs(x, lambda, bw, kerncentres, allownull = TRUE)

  lambda = klambda(bw, "gaussian", lambda)

  z = (x - kerncentres)/lambda
  pnorm(z)
}

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kpuniform <- function(x = 0, lambda = NULL, bw = NULL, kerncentres = 0) {
  
  check.kinputs(x, lambda, bw, kerncentres, allownull = TRUE)
  
  lambda = klambda(bw, "uniform", lambda)
  
  z = pmax(pmin((x - kerncentres)/lambda, 1), -1)
  (z + 1)/2
}

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kptriangular <- function(x = 0, lambda = NULL, bw = NULL, kerncentres = 0) {
  
  check.kinputs(x, lambda, bw, kerncentres, allownull = TRUE)
  
  lambda = klambda(bw, "triangular", lambda)
  
  z = pmax(pmin((x - kerncentres)/lambda, 1), -1)
  ((z + 1)^2 * (z <= 0) + (1 + 2*z - z^2) * (z > 0))/2
}

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kpepanechnikov <- function(x = 0, lambda = NULL, bw = NULL, kerncentres = 0) {
  
  check.kinputs(x, lambda, bw, kerncentres, allownull = TRUE)
  
  lambda = klambda(bw, "epanechnikov", lambda)
  
  z = pmax(pmin((x - kerncentres)/lambda, 1), -1)
  (3*z - z^3 + 2)/4
}

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kpbiweight <- function(x = 0, lambda = NULL, bw = NULL, kerncentres = 0) {
  
  check.kinputs(x, lambda, bw, kerncentres, allownull = TRUE)
  
  lambda = klambda(bw, "biweight", lambda)
  
  z = pmax(pmin((x - kerncentres)/lambda, 1), -1)
  15*(z - 2*z^3/3 + z^5/5 + 8/15)/16
}

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kptriweight <- function(x = 0, lambda = NULL, bw = NULL, kerncentres = 0) {
  
  check.kinputs(x, lambda, bw, kerncentres, allownull = TRUE)
  
  lambda = klambda(bw, "triweight", lambda)
  
  z = pmax(pmin((x - kerncentres)/lambda, 1), -1)
  35*(z - z^3 + 3*z^5/5 - z^7/7 + 16/35)/32
}

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kptricube <- function(x = 0, lambda = NULL, bw = NULL, kerncentres = 0) {
  
  check.kinputs(x, lambda, bw, kerncentres, allownull = TRUE)
  
  lambda = klambda(bw, "tricube", lambda)
  
  z = pmax(pmin((x - kerncentres)/lambda, 1), -1)
  (70*(z + 3*z^4/4 + 3*z^7/7 + z^10/10 + 81/140)/81) * (z <= 0) + 
  (0.5 + 70*(z - 3*z^4/4 + 3*z^7/7 - z^10/10)/81) * (z > 0) 
}


#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kpparzen <- function(x = 0, lambda = NULL, bw = NULL, kerncentres = 0) {
  
  check.kinputs(x, lambda, bw, kerncentres, allownull = TRUE)
  
  lambda = klambda(bw, "parzen", lambda)
  
  z = pmax(pmin((x - kerncentres)/lambda, 1), -1)
  (z < -0.5) * (2 + 8*(z + 1.5*z^2 + z^3 + z^4/4))/3 + 
  ((z >= -0.5) & (z < 0)) * (3 + 8*z - 16*z^3 - 12*z^4)/6 +
  ((z >= 0) & (z <= 0.5)) * (3 + 8*z - 16*z^3 + 12*z^4)/6 + 
  (z > 0.5) * (1 + 8*(z - 1.5*z^2 + z^3 - z^4/4))/3
}

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kpcosine <- function(x = 0, lambda = NULL, bw = NULL, kerncentres = 0) {
  
  check.kinputs(x, lambda, bw, kerncentres, allownull = TRUE)
  
  lambda = klambda(bw, "cosine", lambda)
  
  z = pmax(pmin((x - kerncentres)/lambda, 1), -1)
  (z+1+sin(pi*z)/pi)/2
}

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kpoptcosine <- function(x = 0, lambda = NULL, bw = NULL, kerncentres = 0) {
  
  check.kinputs(x, lambda, bw, kerncentres, allownull = TRUE)
  
  lambda = klambda(bw, "optcosine", lambda)

  z = pmax(pmin((x - kerncentres)/lambda, 1), -1)
  (sin(pi*z/2)+1)/2
}

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kdz <- function(z, kernel = "gaussian") {
  check.kinputs(x = z, lambda = 1, bw = 1, kerncentres = 0)
  check.kernel(kernel)

  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)
  
  switch(kernel,
    gaussian = kdgaussian(z),
    uniform = kduniform(z),
    triangular = kdtriangular(z),
    epanechnikov = kdepanechnikov(z),
    biweight = kdbiweight(z),
    triweight = kdtriweight(z),
    tricube = kdtricube(z),
    parzen = kdparzen(z),
    cosine = kdcosine(z),
    optcosine = kdoptcosine(z))
}

#' @export
#' @aliases kernels kdz kpz
#'  kdgaussian kduniform kdtriangular kdepanechnikov kdbiweight kdtriweight kdtricube kdparzen kdcosine kdoptcosine
#'  kpgaussian kpuniform kptriangular kpepanechnikov kpbiweight kptriweight kptricube kpparzen kpcosine kpoptcosine
#' @rdname kernels
kpz <- function(z, kernel = "gaussian") {
  check.kinputs(x = z, lambda = 1, bw = 1, kerncentres = 0)
  check.kernel(kernel)
  
  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  switch(kernel,
    gaussian = kpgaussian(z),
    uniform = kpuniform(z),
    triangular = kptriangular(z),
    epanechnikov = kpepanechnikov(z),
    biweight = kpbiweight(z),
    triweight = kptriweight(z),
    tricube = kptricube(z),
    parzen = kpparzen(z),
    cosine = kpcosine(z),
    optcosine = kpoptcosine(z))
}
