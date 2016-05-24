#' Local Extrema of Time Series
#' 
#' Find the local minima and maxima from input data. This includes the
#' artificial extrema added to the ends of the data as specified in the
#' original EEMD article [1]. In the case of flat regions at the extrema, the center point of the flat region 
#' will be considered the extremal point [2].
#' 
#' @export
#' @name extrema
#' @param input Numeric vector or time series object.
#' @return a list with matrices \code{minima} and \code{maxima} which give time points and values of local minima and
#' maxima of \code{input} where time points are transformed to match the sampling times of \code{input}.
#' @references
#' \enumerate{
#'  \item{ Z. Wu and N. Huang, "Ensemble Empirical Mode Decomposition: A 
#'       Noise-Assisted Data Analysis Method", Advances in Adaptive Data Analysis,
#'       Vol. 1 (2009) 1--41.}
#'  \item{P. Luukko, J. Helske and E. Räsänen, 
#' "Introducing libeemd: A program package for performing the ensemble empirical mode decomposition", Computational Statistics (2015).}
#' }
#' @examples
#' ext <- extrema(UKgas)
#' plot(UKgas, ylim = range(ext$maxima[, 2], ext$minima[, 2]))
#' points(ext$maxima, col = 2, pch = 19)
#' points(ext$minima, col = 2, pch = 19)
#' 
#' # Artificial extremas obtained by extrapolating last two extrema
#' # Beginning of the series
#' lines(ext$minima[1:3, ], col = 4) 
#' # This is discarded as it produces smaller extrema than the last observation:
#' b <- lm(c(ext$maxima[2:3, 2]) ~ ext$maxima[2:3, 1])$coef[2]
#' points(x = ext$maxima[1, 1], y = ext$maxima[2, 2] - b, col = 4,pch = 19) 
#' lines(x = ext$maxima[1:3, 1], y = c(ext$maxima[2, 2] - b, ext$maxima[2:3, 2]), col = 4)
#' # End of the series
#' # These produce more extreme values than the last observation which is thus disregarded
#' lines(ext$minima[27:29, ],col = 4) 
#' lines(ext$maxima[26:28, ],col = 4) 
#' 
extrema <- function(input) {  
  output <- extremaR(input)
   if (inherits(input, "ts")) {
     output$x_max <- time(input)[output$x_max + 1]
     output$x_min <- time(input)[output$x_min + 1]
   }
  list(minima = cbind(time = output$x_min, value = output$y_min), 
       maxima = cbind(time = output$x_max, value = output$y_max))
}