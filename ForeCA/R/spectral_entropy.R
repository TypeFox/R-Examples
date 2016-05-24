#' @title Estimates spectral entropy of a time series
#' 
#' @description
#' Estimates \emph{spectral entropy} from a univariate (or multivariate) 
#' normalized spectral density.
#' 
#' @details
#' The \emph{spectral entropy} equals the Shannon entropy of the spectral density
#' \eqn{f_x(\lambda)} of a stationary process \eqn{x_t}: 
#' \deqn{ 
#' H_s(x_t) = - \int_{-\pi}^{\pi} f_x(\lambda) \log f_x(\lambda) d \lambda, 
#' }
#' where the density is normalized such that 
#' \eqn{\int_{-\pi}^{\pi} f_x(\lambda) d \lambda = 1}. An estimate of \eqn{f(\lambda)}
#'  can be obtained 
#' by the (smoothed) periodogram (see \code{\link{mvspectrum}}); thus using discrete, and 
#' not continuous entropy.
#' 
#' @inheritParams common-arguments
#' @param series univariate time series of length \eqn{T}.  In the rare case
#' that users want to call this for a multivariate time \code{series}, note 
#' that the estimated spectrum is in general \emph{not} normalized for the computation. 
#' Only if the original data is whitened, then it is normalized.
#' @param mvspectrum.output optional; one can directly provide an estimate of 
#' the spectrum of \code{series}. Usually the output of \code{\link{mvspectrum}}.
#' @param ... additional arguments passed to \code{\link{mvspectrum}}.
#' @return 
#' A non-negative real value for the spectral entropy \eqn{H_s(x_t)}.
#' @seealso \code{\link{Omega}}, \code{\link{discrete_entropy}}
#' @references 
#' Jerry D. Gibson and Jaewoo Jung (2006). \dQuote{The
#' Interpretation of Spectral Entropy Based Upon Rate Distortion Functions}.
#' IEEE International Symposium on Information Theory, pp. 277-281.
#' 
#' L. L. Campbell, \dQuote{Minimum coefficient rate for stationary random processes}, 
#' Information and Control, vol. 3, no. 4, pp. 360 - 371, 1960.
#' 
#' @keywords ts univar math
#' @export
#' @examples
#' 
#' set.seed(1)
#' eps <- rnorm(100)
#' spectral_entropy(eps)
#' 
#' phi.v <- seq(-0.95, 0.95, by = 0.1)
#' kMethods <- c("wosa", "multitaper", "direct", "pgram")
#' SE <- matrix(NA, ncol = length(kMethods), nrow = length(phi.v))
#' for (ii in seq_along(phi.v)) {
#'   xx.tmp <- arima.sim(n = 200, list(ar = phi.v[ii]))
#'   for (mm in seq_along(kMethods)) {
#'     SE[ii, mm] <- spectral_entropy(xx.tmp, spectrum.control = 
#'                                      list(method = kMethods[mm]))
#'   }
#' }
#' 
#' matplot(phi.v, SE, type = "l", col = seq_along(kMethods))
#' legend("bottom", kMethods, lty = seq_along(kMethods), 
#'        col = seq_along(kMethods))
#'        
#' # AR vs MA
#' SE.arma <- matrix(NA, ncol = 2, nrow = length(phi.v))
#' SE.arma[, 1] <- SE[, 2]
#' 
#' for (ii in seq_along(phi.v)){
#'   yy.temp <- arima.sim(n = 1000, list(ma = phi.v[ii]))
#'   SE.arma[ii, 2] <- 
#'     spectral_entropy(yy.temp, spectrum.control = list(method = "multitaper"))
#' }
#' 
#' matplot(phi.v, SE.arma, type = "l", col = 1:2, xlab = "parameter (phi or theta)",
#'         ylab = "Spectral entropy")
#' abline(v = 0, col = "blue", lty = 3)
#' legend("bottom", c("AR(1)", "MA(1)"), lty = 1:2, col = 1:2)
#' 

spectral_entropy <- function(series = NULL, spectrum.control = list(),
                             entropy.control = list(),
                             mvspectrum.output = NULL, ...){
  
  stopifnot(xor(is.null(series), is.null(mvspectrum.output)),
            is.null(mvspectrum.output) || 
              length(dim(mvspectrum.output)) >= 0)

  if (is.null(mvspectrum.output)) {
    series <- as.matrix(series)
    is.whitened <- (isTRUE(all.equal(cov(series), diag(1, ncol(series)))) && 
                      isTRUE(all.equal(colMeans(series), 0)))
    attr(series, "whitened") <- is.whitened
    num.series <- ncol(series)
    spectrum.control <- complete_spectrum_control(spectrum.control)
    mvspectrum.output <- mvspectrum(series, 
                                    method = spectrum.control$method, 
                                    smoothing = spectrum.control$smoothing,
                                    normalize = attr(series, "whitened"),
                                    ...)
  }
  if (is.null(dim(mvspectrum.output))) {
    dim(mvspectrum.output) <- c(length(mvspectrum.output), 1, 1)
  }
  num.series <- dim(mvspectrum.output)[2]
  
  if (num.series > 1) {
    num.freqs <- dim(mvspectrum.output)[1]
  } else {
    num.freqs <- length(mvspectrum.output)
  }
  
  entropy.control$base <- NULL
  entropy.control <- complete_entropy_control(entropy.control,
                                              num.outcomes = num.freqs * 2)
  if (num.series > 1) {
    spec.ent <- mvspectrum2wcov(mvspectrum.output * 
                                  -log(mvspectrum.output, 
                                       base = entropy.control$base))
  } else {
    # this has sum() = 1
    spec.dens <- c(rev(c(mvspectrum.output)), c(mvspectrum.output))
    spec.dens[spec.dens < 0] <- 0
    spec.dens <- spec.dens / sum(spec.dens)
    spec.ent <- do.call("discrete_entropy", 
                        c(list(probs = spec.dens), entropy.control))
  }
  attr(spec.ent, "nfrequencies") <- num.freqs
  attr(spec.ent, "base") <- entropy.control$base
  return(spec.ent)
} 