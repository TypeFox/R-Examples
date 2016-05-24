#' @title S3 methods for class 'mvspectrum'
#' @name mvspectrum-utils
#' @description
#' S3 methods for multivariate spectrum estimation. 
#' @param x an object of class \code{"foreca.one_weightvector"}.
#' @param ... additional arguments passed to \code{\link[graphics]{matplot}}.
#' @seealso \code{\link{get_spectrum_from_mvspectrum}}
#' @examples
#' # see examples in 'mvspectrum'
#'
NULL


#' @rdname mvspectrum-utils
#' @method plot mvspectrum
#' @description
#' \code{plot.mvspectrum} plots all univariate spectra. Analogouos to 
#' \code{\link[stats]{spectrum}} when \code{plot = TRUE}.
#' @keywords manip hplot
#' @param log logical; if \code{TRUE} (default), it plots the spectra on 
#' log-scale.
#' @export
#' 
#' @examples
#' SS <- mvspectrum(diff(log(EuStockMarkets)) * 100, 
#'                  spectrum.control = list(method = "multitaper"))
#' plot(SS, log = FALSE)
#' 

plot.mvspectrum <- function(x, log = TRUE, ...) {
  
  freqs <- attr(x, "frequency")
  uni.spectra <- get_spectrum_from_mvspectrum(x)
  
  matplot(freqs, uni.spectra, type = "l", lty = 1, 
          xlab = bquote(lambda), ylab = bquote(hat(f)(lambda)), 
          log = ifelse(log, "y", ""), ...)
} 

