#' Mean Residual Life Plot for the Generalized Pareto Distribution
#'
#' Plots the empirical mean residual life, with confidence intervals. The mean residual life plot provides a
#' visual diagnostic tool to choose a threshold for exceedances.
#' @param data Vector of data.
#' @param thresholds A numeric vector of threshold(s) to plot vertically. Defaults to NULL.
#' @param conf The level of the confidence bounds to use. Defaults to 0.95.
#' @param npoints The number of points to interpolate with. Defaults to 100.
#' @examples
#' ## Not run
#' ## x <- rgpd(500, loc = 0, scale = 1, shape = 0.1)
#' ## mrlPlot(x, thresholds = c(2))
#' @export

mrlPlot <- function(data, thresholds = NULL, conf = .95, npoints = 100) {
  umin <- min(data)
  umax <- findthresh(data, 2)
  y <- yu <- yl <- rep(0, npoints)
  u <- seq(umin, umax, length = npoints)
  for(i in 1:npoints) {
    excess <- data[data > u[i]]
    y[i] <- mean(excess - u[i])
    sdev <- sqrt(var(excess))
    n <- length(excess)
    yu[i] <- y[i] + (qnorm((1 + conf)/2) * sdev)/sqrt(n)
    yl[i] <- y[i] - (qnorm((1 + conf)/2) * sdev)/sqrt(n)
  }
  plot(u, y, type = "l", xlab = "Threshold", ylab = "Mean Excess",
       ylim = c(min(yl), max(yu)))
  lines(u, yl, lty = 2, col = 4)
  lines(u, yu, lty = 2, col = 4)
  if(!is.null(thresholds)) {
    for(j in 1:length(thresholds)) {
      abline(v = thresholds[j], untf = FALSE, col = j)
    }
  }
  points(data, rep(min(yl), length(data)))
}


