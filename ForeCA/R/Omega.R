#' @title Estimate forecastability of a time series
#' 
#' @description
#' An estimator for the forecastability \eqn{\Omega(x_t)} of a univariate time series \eqn{x_t}.
#' Currently it uses a discrete plug-in estimator given the empirical spectrum (periodogram).
#' 
#' @details
#' The \emph{forecastability} of a stationary process \eqn{x_t} is defined as 
#' (see References)
#' 
#' \deqn{ 
#' \Omega(x_t) = 1 - \frac{ - \int_{-\pi}^{\pi} f_x(\lambda) \log f_x(\lambda) d \lambda }{\log 2 \pi} \in [0, 1]
#' }
#' where \eqn{f_x(\lambda)} is the normalized spectral \emph{density} of \eqn{x_t}. 
#' In particular \eqn{ \int_{-\pi}^{\pi} f_x(\lambda) d\lambda = 1}.
#'
#' 
#' For white noise \eqn{\varepsilon_t} forecastability 
#' \eqn{\Omega(\varepsilon_t) = 0}; for a sum of sinusoids it equals \eqn{100} \%. 
#' However, empirically it reaches \eqn{100\%} only if the estimated spectrum has
#' exactly one peak at some \eqn{\omega_j} and \eqn{\widehat{f}(\omega_k) = 0}
#'  for all \eqn{k\neq j}.
#' 
#' In practice, a time series of length \code{T} has \eqn{T} Fourier frequencies
#' which represent a discrete 
#' probability distribution.  Hence entropy of \eqn{f_x(\lambda)} must be 
#' normalized by \eqn{\log T}, not by \eqn{\log 2 \pi}.
#' 
#' Also we can use several smoothing techniques to obtain a less variance estimate of 
#' \eqn{f_x(\lambda)}.
#' 
#' @param series a univariate time series; if it is multivariate, then 
#' \code{\link{Omega}} works component-wise (i.e., same as \code{apply(series, 2, Omega)}).
#' @inheritParams common-arguments
#' @inheritParams spectral_entropy
#' @return 
#' A real-value between \eqn{0} and \eqn{100} (\%). \eqn{0} means not
#' forecastable (white noise); \eqn{100} means perfectly forecastable (a
#' sinusoid).
#' @export
#' @seealso \code{\link{spectral_entropy}}, \code{\link{discrete_entropy}}, 
#' \code{\link{continuous_entropy}}
#' @references Goerg, G. M. (2013). \dQuote{Forecastable Component
#'  Analysis}. Journal of Machine Learning Research (JMLR) W&CP 28 (2): 64-72, 2013.
#'  Available at \url{jmlr.org/proceedings/papers/v28/goerg13.html}.
#' @keywords math univar
#' @examples
#' 
#' nn <- 100
#' eps <- rnorm(nn)  # white noise has Omega() = 0 in theory
#' Omega(eps, spectrum.control = list(method = "direct")) 
#' # smoothing makes it closer to 0
#' Omega(eps, spectrum.control = list(method = "wosa"))
#' 
#' xx <- sin(seq_len(nn) * pi / 10)
#' Omega(xx, spectrum.control = list(method = "direct"))
#' Omega(xx, entropy.control = list(threshold = 1/40)) 
#' Omega(xx, spectrum.control = list(method = "wosa"), 
#'       entropy.control = list(threshold = 1/20))
#' 
#' # an AR(1) with phi = 0.5
#' yy <- arima.sim(n = nn, model = list(ar = 0.5))
#' Omega(yy, spectrum.control = list(method = "wosa"))
#' 
#' # an AR(1) with phi = 0.9 is more forecastable
#' yy <- arima.sim(n = nn, model = list(ar = 0.9))
#' Omega(yy, spectrum.control = list(method = "wosa"))
#' 
Omega <- function(series = NULL, 
                  spectrum.control = list(),
                  entropy.control = list(),
                  mvspectrum.output = NULL) {

  stopifnot(xor(is.null(series), is.null(mvspectrum.output)))
  
  if (is.null(mvspectrum.output)) {
    series <- as.matrix(series)
    # center and scale the data
    series <- scale(series)
    num.series <- ncol(series)
    num.freqs <- floor(nrow(series) / 2)
  } else {
    num.freqs <- nrow(as.matrix(mvspectrum.output))
    num.series <- ncol(as.matrix(mvspectrum.output))
  }
  
  entropy.control$base <- NULL
  entropy.control <- complete_entropy_control(entropy.control,
                                              num.outcomes = 2 * num.freqs)
  spectrum.control <- complete_spectrum_control(spectrum.control)
  
  if (num.series > 1) {
    stopifnot(is.null(mvspectrum.output))
    OMEGAs <- apply(series, 2, 
                    function(x) {
                      attr(x, "whitened") <- attr(series, "whitened")
                      Omega(x, spectrum.control = spectrum.control, 
                            entropy.control = entropy.control, 
                            mvspectrum.output = mvspectrum.output)
                    })
    attr(OMEGAs, "unit") <- "%"
    return(OMEGAs)
  } else {
    h.spectral <- spectral_entropy(series = series, 
                                   spectrum.control = spectrum.control,
                                   mvspectrum.output = mvspectrum.output, 
                                   entropy.control = entropy.control)
    #print(h.spectral)
    OMEGA <- (1 - c(h.spectral)) * 100
    attr(OMEGA, "unit") <- "%"
    return(OMEGA)
  }
} 