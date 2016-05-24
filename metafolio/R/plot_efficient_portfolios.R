#' Basic plot of efficient portfolio and asset contributions
#'
#' This function creates a mean-variance plot of the portfolios across possible
#' asset weights, colour the efficient frontier, and show the contribution of
#' the different stocks/assets. It also (invisibly) returns the values that make
#' up the plot so you can create your own custom plots with the data. See the
#' Returns section for more details.
#'
#' @param port_vals A matrix of means and variances (down the two columns).
#'   This likely comes from the output of \code{\link{monte_carlo_portfolios}}.
#' @param weights_matrix The same weight matrix that was passed to
#'   \code{\link{monte_carlo_portfolios}}.
#' @param pal Colour palette for the stocks/assets in the barplot.
#' @param plot Logical: should the plots be made?
#' @param ylab_dots Y axis label for the mean-variance scatterplot.
#' @param xlab_dots X axis label for the mean-variance scatterplot.
#' @param ylab_bars Y axis label for the barplot.
#' @param xlab_bars X axis label for the barplot.
#' @param port_cols Colours for the dots. A vector of colours for the
#'   non-efficient and efficient portfolios.
#' @param pch Dot type
#' @param ... Anything else to pass to both
#'   \code{\link[graphics]{plot.default}} and \code{\link[graphics]{barplot}}.
#' @return
#' A two panel plot and an (invisible) list of values calculated within the
#' function. This list contains \code{pv} (mean, variance, and whether it was
#' part of the efficient frontier); \code{ef_port_ids} (the portfolio IDs [run
#' numbers] that are part of the efficient frontier; \code{min_var_port_id} (the
#' portfolio ID for the minimum-variance portfolio); \code{ef_weights} (the
#' weights of the portfolios on the efficient frontier).
#' @export
#' @examples
#' \dontrun{
#' weights_matrix <- create_asset_weights(n_pop = 6, n_sims = 3000,
#' weight_lower_limit = 0.001)
#' mc_ports <- monte_carlo_portfolios(weights_matrix = weights_matrix,
#'  n_sims = 3000, mean_b = 1000)
#' 
#' col_pal <- rev(gg_color_hue(6))
#' ef_dat <- plot_efficient_portfolios(port_vals = mc_ports$port_vals,
#'  pal = col_pal, weights_matrix = weights_matrix)
#' names(ef_dat)
#' }

plot_efficient_portfolios <- function(port_vals, weights_matrix, pal,
  plot = TRUE, ylab_dots = "Mean of metapopulation growth rate",
  xlab_dots = "Variance of metapopulation growth rate", ylab_bars =
  "Percentage", xlab_bars = "Variance (multiplied by 1000)", port_cols
  = c("grey50", "red"), pch = 19, ...) {

  pv <- as.data.frame(port_vals)
  names(pv) <- c("m", "v")
  # let's pull out the efficient frontier:
  m.bins <- seq(min(pv$m), max(pv$m), length.out = 50)
  pv$m.bin <- findInterval(pv$m, m.bins)
  pv$optim_set <- 0
  for(i in unique(pv$m.bin)) {
    pv_to_check <- which(pv$m.bin == i)
    pv[pv_to_check, "optim_set"][which.min(pv[pv_to_check, ]$v)] <- 1
  }

  # those below the min variance port should be ignored:
  # they are not desirable
  m_at_min_var_port <- pv$m[pv$v == min(pv$v)]
  pv$optim_set[pv$m < m_at_min_var_port] <- 0

  pv$id <- 1:nrow(pv)
  pv <- pv[order(pv$optim_set), ]

  if(plot == TRUE){
    par(mfrow = c(1, 2), xpd = NA)
    par(cex = 0.8)
    with(pv, plot(v, m, pch = pch, col = port_cols[optim_set+1], cex =
        0.8, xlab = xlab_dots, ylab = ylab_dots, ...))
  }

  pv <- pv[order(pv$v), ]
  # efficient frontier portfolios
  ef_ports <- pv[pv$optim_set == 1, "id"]
  ef_weights <- weights_matrix[ef_ports, ]

  if(plot == TRUE){
    barplot(t(100*ef_weights), col = pal,border = "grey50", names.arg
      = round(pv$v[pv$id %in% ef_ports]*1000, 2), xlab = xlab_bars,
      ylab = ylab_bars, las =1, ...)
  }
  invisible(list(pv = pv, ef_port_ids = ef_ports, min_var_port_id =
      ef_ports[1], ef_weights = ef_weights))
}
