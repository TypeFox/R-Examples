#' @rdname LSC-utils
#' @description
#' \code{\link{plot_LSC_0plus1D}} plots LSC for a (0+1)D field, i.e., a 
#' time series.
#' @keywords hplot
#' @export
#' @examples
#' 
#' state.sim <- rpois(100, 1)
#' 
#' lsc.est <- states2LSC(states = state.sim)
#' class(lsc.est) <- c("LSC", "LSC_0plus1D")
#' plot_LSC_0plus1D(lsc.est)
#' 
#' weights.sim <- matrix(runif(1000, 0, 1), ncol = 10)
#' weights.sim <- normalize(weights.sim)
#' lsc.est <- states2LSC(weight.matrix = weights.sim)
#' plot_LSC_0plus1D(lsc.est)
#' 

plot_LSC_0plus1D = function(z, col = NULL, lsc.unit = "bits", ...) {
  
  data.tmp <- data.frame(lsc = z[],
                         time = seq_along(z[]))
  data.tmp$smoothed <- loess(lsc ~ time, span = 0.25, data = data.tmp)$fit
  
  matplot(data.tmp$time, 
          data.tmp[, c("lsc", "smoothed")], type = "l",
          lty = c(1, 2), col = c(1, 2), lwd = c(1, 2),
          ylab = lsc.unit, xlab = 'Time', ...)
}

