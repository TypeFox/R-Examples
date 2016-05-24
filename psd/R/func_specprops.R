#' Calculate properties of multitaper power spectral density estimates
#'
#' @description
#' Various spectral properties may be computed from the vector of tapers, and
#' if necessary the sampling frequency.
#' 
#' @details
#' Parameter Details:
#' \subsection{Uncertainty}{
#' See \code{\link{spec_confint}} for details.
#' }
#'
#' \subsection{Resolution}{
#' The frequency resolution depends on the number of tapers (\eqn{K}), and
#' is found from 
#' \deqn{\frac{K \cdot f_N}{N_f}} 
#' where \eqn{f_N} is the Nyquist
#' frequency and \eqn{N_f} is the 
#' number of frequencies estimated.
#' }
#'
#' \subsection{Degrees of Freedom}{
#' There are two degrees of freedom for each taper \eqn{K}:
#' \deqn{\nu = 2 K}
#' }
#'
#' \subsection{Bandwidth}{
#' The bandwidth of a multitaper estimate depends on the number of
#' tapers.
#' Following Walden et al (1995) the effective bandwidth is \eqn{\approx 2W}
#' where
#' \deqn{W = \frac{K + 1}{2 N}} 
#(N+1)}}
#' and \eqn{N} is the number of terms in the series, which makes \eqn{N \cdot W} the
#' approximate time-bandwidth product.
#' }
#' 
#' @author A.J. Barbour
#' @name spectral_properties
#' @export
#' @seealso \code{\link{spec_confint}}, \code{\link{psd-package}}
#'
#' @param tapvec object with class \code{'tapers'} or \code{'spec'}
#' @param f.samp scalar; the sampling frequency (e.g. Hz) of the series the tapers are for
#' @param n.freq scalar; the number of frequencies of the original spectrum 
#' (if \code{NULL} the length of the tapers object is assumed to be the number)
#' @param p numeric; the coverage probability, bound within \eqn{[0,1)}
#' @param  db.ci logical; should the uncertainty confidence intervals be returned as decibels?
#' @param ... additional arguments (unused)
#' @return A list with the following properties (and names):
#' \itemize{
#' \item{\code{taper}: The original taper vector.}
#' \item{\code{stderr.chi .upper, .lower, .median}: results returned from \code{\link{spec_confint}}.}
#' \item{\code{resolution}: The effective spectral resolution.}
#' \item{\code{dof}: The number of degrees of freedom.}
#' \item{\code{bw}: The effective bandwidth of the spectrum.}
#' }
#'
#' @example inst/Examples/rdex_spectralproperties.R
spectral_properties <- function(tapvec, ...) UseMethod("spectral_properties")

#' @rdname spectral_properties
#' @aliases spectral_properties.spec
#' @export
spectral_properties.spec <- function(tapvec, ...){
  stopifnot(is.spec(Pspec <- tapvec))
  n.freq <- length(Pspec[['freq']])
  f.samp <- 2 * Pspec[['freq']][n.freq]
  tapvec <- Pspec[['taper']]
  spectral_properties(tapvec, f.samp, n.freq, ...)
}

#' @rdname spectral_properties
#' @aliases spectral_properties.tapers
#' @export
spectral_properties.tapers <- function(tapvec, ...){
  stopifnot(is.tapers(tapvec))
  K <- as.integer(tapvec)
  spectral_properties(K, ...)
}

#' @rdname spectral_properties
#' @aliases spectral_properties.default
#' @export
spectral_properties.default <- function(tapvec, f.samp=1, n.freq=NULL, p=0.95, db.ci=FALSE, ...){
  
  K <- as.vector(tapvec)
  if (is.null(n.freq)) n.freq <- length(K)
  
  Nyquist <- f.samp/2
  
  #Var <- 10 / K / 12
  ## Deg Freedom: PW93 Ch7 343
  Dof <- 2 * K
  #
  ##BW <- K / n.freq
  ##Resolu <- BW * Nyquist
  #
  ## Bandwidth
  # Walden et al
  # half-width W = (K + 1)/{2(N + 1)}
  # effective bandwidth ~ 2 W (accurate for many spectral windows)
  W <- (K+1)/(2*n.freq)
  BW <- 2 * W
  ## Resolution
  Resolu <- 2 * BW
  ## Uncertainty CI -- how does it compare to StdErr
  StdErrCI <- spec_confint(Dof, p, as.db=db.ci)
  ##
  return(data.frame(taper=K, stderr.chi=StdErrCI, resolution=Resolu, dof=Dof, bw=BW))
}

#' Confidence intervals for multitaper power spectral density estimates 
#'
#' @details
#' The errors are estimated 
#' from the number of degrees of freedom \eqn{\nu} by evaluating
#' the \eqn{\chi_{p,\nu}^{2}(\nu,\nu)} distribution for an optional 
#' coverage probability \eqn{p} (defaulting to \eqn{p=0.95}).  
#' Additionally, the
#' \eqn{p=0.5} values and an approximation from \eqn{1/\sqrt{\nu - 1}}
#' are returned.
#'
#' A more 
#' sophisticated (and complicated) approach would be to
#' estimate via jack-knifing (Prieto et al 2007), but this is not yet
#' made available.
#'
#' Additive uncertainties \eqn{\delta S} are returned, such that 
#' the spectrum with confidence interval is \eqn{S \pm \delta S}.
#'
#' @author A.J. Barbour; some code modified from the \code{spec.ci} function inside \code{plot.spec}
#' @name spec_confint
#' @export
#' @seealso \code{\link{spectral_properties}}, \code{\link{psd-package}}, \code{plot.spec}, \code{\link{dB}}
#' @param dof numeric; the degrees of freedom \eqn{\nu}
#' @param p numeric; the coverage probability \eqn{p}, bound within \eqn{[0,1)}
#' @param as.db logical; should the values be returned as decibels?
#' @return A \code{data.frame} with the following properties (and names):
#' \itemize{
#' \item{\code{lower}: Based on upper tail probabilities (\eqn{p})}
#' \item{\code{upper}: Based on lower tail probabilities (\eqn{1-p})}
#' \item{\code{median}: Based on lower tail probabilities (\eqn{p=0.5})}
#' \item{\code{approx}: Approximation based on \eqn{1/\sqrt(\nu - 1)}.}
#' }
#' @example inst/Examples/rdex_confint.R
spec_confint <- function(dof, p = 0.95, as.db=FALSE) UseMethod("spec_confint")

#' @rdname spec_confint
#' @aliases spec_confint.spec
#' @export
spec_confint.spec <- function(dof, p = 0.95, as.db=FALSE){
  stopifnot(is.spec(dof))
  dof <- dof$df
  spec_confint(dof, p, as.db)
}

#' @rdname spec_confint
#' @aliases spec_confint.tapers
#' @export
spec_confint.tapers <- function(dof, p = 0.95, as.db=FALSE){
  stopifnot(is.tapers(dof))
  # two degrees of freedom per taper 
  dof <- 2 * unclass(dof)
  spec_confint(dof, p, as.db)
}

#' @rdname spec_confint
#' @aliases spec_confint.default
#' @export
spec_confint.default <- function(dof, p = 0.95, as.db=FALSE) {
  # Mostly from spec.ci, lifted from plot.spec
  if (p < 0 || p >= 1) stop("coverage probability out of range [0,1)")
  ptail <- (1 - p)
  # qchisq gives distribution function
  # if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
  upper.qp <- 1 - ptail * pchisq(dof, dof, lower.tail = FALSE)
  med.qp <-  0.5 * pchisq(dof, dof, lower.tail = TRUE)
  lower.qp <- ptail * pchisq(dof, dof, lower.tail = TRUE)
  ndof <- length(dof)
  # qchisq gives quantile function
  # spec.ci calculates with 1/(spec/dof)
  ci.ul <- 1 / (qchisq(c(upper.qp, lower.qp, med.qp), dof)/dof)
  # dS = S*ci.ul/dof
  ## heuristically tuned to approximate the median distribution of Chi^2 uncertainties
  approx <- 1 + 1/sqrt(dof-1)
  ci <- data.frame(lower=ci.ul[1:ndof], 
                   upper=ci.ul[(ndof+1):(2*ndof)], 
                   median=ci.ul[(2*ndof+1):(3*ndof)],
                   approx=approx)
  if (as.db) ci <- dB(ci)
  return(ci)
}
