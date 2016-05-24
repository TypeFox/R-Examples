gpdRlPlot <- function(z, conf = 0.95, method = c("delta", "profile")) {
  if(!z$stationary)
    stop("Model must be stationary")
  method <- match.arg(method)
  emp <- z$n / z$npp
  p <- c(seq(from = 0.001, to = emp, by = emp / 100), seq(from = emp, to = 100 * emp, by = emp))
  levels <- matrix(0, length(p), 3)
  for(i in 1:nrow(levels)) {
    y <- gpdRl(z, p[i], conf = conf, method = method)
    levels[i, 1] <- y$Estimate
    levels[i, 2:3] <- y$CI
  }
  plot((z$n + 1) / ((1:z$n) * z$npp), rev(sort(z$data)), type = "n", log = "x",
       xlab = "Return Period", ylab = "Return Level", xlim = c(min((z$n + 1) / ((1:z$n) * z$npp)), max(p)),
       ylim = c(min(z$data, levels[, 2]), max(z$data, levels[, 3])))
  title("Return Level Plot")
  lines(p, levels[, 1])
  lines(p, levels[, 2], col = 4)
  lines(p, levels[, 3], col = 4)
  points((z$n + 1) / ((1:z$n) * z$npp), rev(sort(z$data)))
}


gpdHist <- function(z) {
  if(!z$stationary)
    stop("Model must be stationary")
  excess <- z$data[z$data > z$threshold]
  h <- hist(excess, plot = FALSE)
  x <- seq(min(excess), max(excess), (max(excess) - min(excess))/1000)
  y <- dgpd(x, loc = z$threshold, scale = z$par.ests[1], shape = z$par.ests[2])
  hist(excess, freq = FALSE, xlim = c(min(excess), max(excess)), ylim = c(0, max(max(h$density), max(y))),
       xlab = "x", ylab = "Density", main = "Density Plot")
  points(excess, rep(0, length(excess)))
  lines(x, y, col = 4)
}


gpdPP <- function(z, scalevec, shapevec) {
  excess <- (z$data - z$threshold)[z$data > z$threshold]
  n <- length(excess)
  Series <- seq(1, n, 1) / (n+1)
  p <- pgpd(excess, loc = 0, scale = scalevec, shape = shapevec)
  plot(sort(p), Series, xlab = "Empirical", ylab = "Model", xlim = c(0,1), ylim = c(0,1))
  if(z$stationary)
    title("Probability Plot")
  if(!z$stationary)
    title("Residual Probability Plot")
  abline(0, 1, col = 4)
}


gpdQQ <- function(z, scalevec, shapevec) {
  excess <- (z$data - z$threshold)[z$data > z$threshold]
  n <- length(excess)
  Series <- seq(1, n, 1) / (n+1)
  emp <- qgpd(Series, loc = 0, scale = scalevec, shape = shapevec)
  plot(sort(excess), emp, xlab = "Empirical", ylab = "Model",
       xlim = c(min(excess, emp), max(excess, emp)), ylim = c(min(excess, emp), max(excess, emp)))
  if(z$stationary)
    title("Quantile Plot")
  if(!z$stationary)
    title("Residual Quantile Plot")
  abline(0, 1, col = 4)
}


## Plots residuals vs. the covariates
gpdResid <- function(z, scalevec, shapevec) {
  if(z$stationary)
    stop("Model cannot be stationary")
  excess <- (z$data - z$threshold)[z$data > z$threshold]
  resid <- pgpd(excess, loc = 0, scale = scalevec, shape = shapevec)
  if(z$parnum[1] > 1) {
    for(i in 2:z$parnum[1]) {
      plot(z$covars[[1]][, i], resid, xlab = paste("Scale", colnames(z$covars[[1]])[i], sep = " "), ylab = "Residuals")
      lines(lowess(z$covars[[1]][, i], resid), col = "red")
    }
  }
  if(z$parnum[2] > 1) {
    for(i in 2:z$parnum[2]) {
      plot(z$covars[[2]][, i], resid, xlab = paste("Shape", colnames(z$covars[[2]])[i], sep = " "), ylab = "Residuals")
      lines(lowess(z$covars[[2]][, i], resid), col = "red")
    }
  }
}


#' Diagnostic plots for a fit to the Generalized Pareto distribution
#'
#' @param z A class object returned from `gpdFit'.
#' @param conf Confidence level used in the return level plot.
#' @param method The method to compute the return level confidence interval - either delta method (default) or
#' profile likelihood. Choosing profile likelihood may be quite slow.
#' @examples
#' ## Not run
#' # x <- rgpd(10000, loc = 0.5, scale = 1, shape = 0.1)
#' # z <- gpdFit(x, nextremes = 500)
#' # plot(z)
#' @return For stationary models, provides return level, density, probability, and quantile plots for the GPD exceedances. The
#' overlaid density is the `true' density for the estimated parameters. For nonstationary models, provides
#' residual probability and quantile plots. In addition, nonstationary models provide plots of the residuals vs.
#' the parameter covariates.
#' @details See the reference for details on how return levels are calculated.
#' @references Coles, S. (2001). An introduction to statistical modeling of extreme values (Vol. 208). London: Springer.
#' @importFrom utils menu
#' @export
gpdDiag <- function(z, conf = 0.95, method = c("delta", "profile")) {
  method <- match.arg(method)
  par(ask = TRUE, mfcol = c(2, 2))
  scalevec <- z$links[[1]](rowSums(t(z$par.ests[1:z$parnum[1]] * t(z$covars[[1]]))))
  shapevec <- z$links[[2]](rowSums(t(z$par.ests[(z$parnum[1] + 1):(z$parnum[1] + z$parnum[2])] * t(z$covars[[2]]))))
  choice <- 1
  while(choice > 0) {
    choice <- menu(c("Return Level Plot", "Density Plot", "PP Plot", "QQ Plot",
                     "Residual Scatterplot(s)"), title = "\nMake a plot selection (or 0 to exit):")
    switch(choice + 1,
           cat("Exited\n"),
           if(!z$stationary) stop("Model must be stationary") else try(gpdRlPlot(z, conf, method), silent = TRUE),
           if(!z$stationary) stop("Model must be stationary") else try(gpdHist(z), silent = TRUE),
           try(gpdPP(z, scalevec, shapevec), silent = TRUE),
           try(gpdQQ(z, scalevec, shapevec), silent = TRUE),
           if(z$stationary) stop("Model cannot be stationary") else try(gpdResid(z, scalevec, shapevec), silent = TRUE)
    )
  }
  par(mfrow = c(1, 1))
}










