#' Plot correlation of returns (i.e. metapopulation abundance) across stocks.
#'
#' Create a matrix plot showing the correlation between the log returns of each
#' stock/asset.
#'
#' @param x A list output object from \code{\link{meta_sim}}.
#' @param burn Number of years to discard at start as burn in.
#' @param pal Colours to label each stock/asset.
#' @param xlab X axis label
#' @param ylab Y axis label
#' @examples
#' arma_env_params <- list(mean_value = 16, ar = 0.1, sigma_env = 2, ma = 0)
#' base1 <- meta_sim(n_pop = 10, env_params = arma_env_params, env_type =
#'   "arma", assess_freq = 5)
#' plot_correlation_between_returns(base1)
#' @export

plot_correlation_between_returns <- function(x, burn = 1:30, pal =
  rev(gg_color_hue(x$n_pop)), xlab = "log of return abundance by population",
  ylab = "log of return abundance by population") {

par(mfrow = c(10, 10), cex = 0.5, mar = c(0,0,0,0), oma = c(4, 4, 0, 0))
for (i in 1:10) {
  for (j in 1:10) {
    if (j < i){ plot(x$A[-burn, i], x$A[-burn, j], log = "xy", axes =
      FALSE, pch = 20, col = "grey20")
    box(col = "grey50")
    }else{
      plot(1,1, type = "n", xlab = "", ylab = "", axes = FALSE)
    }
    if (i == 10) mtext(j, side = 1, line = 1, cex = 0.8, col = pal[j])
    if (j == 1) mtext(i, side = 2, las = 1, line = 1, cex = 0.8, col = pal[i])
  }
}
    mtext(xlab, side = 1, line = 2.5, cex = 0.8, outer = TRUE)
    mtext(ylab, side = 2, line = 2.5, cex = 0.8, outer = TRUE, las = 0)

}
