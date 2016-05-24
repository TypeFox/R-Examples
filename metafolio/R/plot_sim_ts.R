#' Add annotations to panel
#'
#' @param label The text to add as a label
#' @param xfrac Fraction over from the left
#' @param yfrac Fraction down from the top
#' @param pos Position of text to pass to \code{\link[graphics]{text}}
#' @param cex Character expansion value to pass to \code{\link[graphics]{text}}
#' @param ... Anything else to pass to \code{\link[graphics]{text}}
annotate <- function(label, xfrac = 0.008, yfrac = 0.18, pos = 4, cex = 0.9, ...) {
  u <- par("usr")
  x <- u[1] + xfrac * (u[2] - u[1])
  y <- u[4] - yfrac * (u[4] - u[3])
  text(x, y, label, pos = pos, cex = cex, ...)
}

#' Add a pretty axis
#'
#' @param side Number indicating the side to add an axis (as in the side number
#' passed to \code{\link[graphics]{axis}}).
#' @param shade_years An optional numerical vector of length two giving the
#' minimum and maximum years over which to add a light grey shading.
#' @param ylab Y axis label
#' @param yticks Logical: should y-axis ticks be added?
my.axis <- function(side, shade_years = NULL, ylab = "", yticks = NA) {
  if(!is.null(shade_years)) {
    rect(min(shade_years), -100, max(shade_years), 1e9, col =
      "#00000030", border = NA)
  }
  #axis(side, col = "grey50")
  if(is.na(yticks[1])) {
    axis(side, col = "grey50", at = pretty(axTicks(2), n= 2))
  } else {
    axis(side, col = "grey50", at = yticks)
  }
  par(las = 0)
  mtext(ylab, side = 2, cex = 0.7, line = 3)
  par(las = 1)
}


#' Standard matrix plot of values by stream for one panel:
#'
#' @param dat The matrix of values to plot
#' @param ymin Minimum y value for axis
#' @param ystretch A fraction to multiply the max value of when setting the y
#' axis limits. This is useful to make space for a panel label within the plot.
#' @param ... Anything else to pass to \code{\link[graphics]{matplot}}.
plot_panel_lines <- function(dat, ymin = c("zero", "min"), ystretch = 1.1, ...) {
  dat[is.na(dat)] <- 0
  if(ymin[1] == "zero")
    ylim <- c(0, max(dat) * ystretch)
  if(ymin[1] == "min")
    ylim <- c(min(dat), max(dat) * ystretch)
  matplot(dat, type = "l", lty = 1, xlab =
    "Time", xaxt = "n", axes = FALSE, xaxs = "i", ylim = ylim , ...)
  box(col = "grey50")
}


#' Plot various time series from a simulation run
#'
#' This function lets you quickly visualize the time series of output from a
#' simulation run.
#'
#' @param x A list output object from a simulation run of
#'   \code{link{meta_sim}}.
#' @param pal A colour palette for the lines. One colour per line (each
#' line is a population time series).
#' @param years_to_show How many years to plot after the burn in period.
#' @param burn The number of years to discard as burn in at the beginning of
#'   the time series.
#' @param adj \code{adj} parameter to pass to \code{\link[graphics]{mtext}} for
#'   panel labels
#' @param shade_years Shade some years? Give a vector. Shading will be applied
#' from the minimum to maximum value. Can be used to show burn in period.
#' @param add_units Should the units be added to the y axis?
#' @param yticks Position of ticks on the Y axis.
#' @param oma \code{oma} vector to pass to \code{par} for outer margin space.
#' @export
#' @examples
#' arma_env_params <- list(mean_value = 16, ar = 0.1, sigma_env = 2, ma = 0)
#' base1 <- meta_sim(n_pop = 10, env_params = arma_env_params, env_type =
#'   "arma", assess_freq = 5, decrease_b = 10)
#' plot_sim_ts(base1, years_to_show = 70, burn = 1:30)

plot_sim_ts <- function(x, pal = rev(gg_color_hue(x$n_pop)),
  years_to_show = 30, burn = 1:50, shade_years = NULL, adj = 0.02,
  add_units = FALSE, yticks = rep(list(NA), 10), oma = c(4, 4.5, 1, 1)) {

  if(!add_units) {
    ylabs <- rep("", 10)
  } else {
    ylabs <- c(rep("Value", 2), rep("#", 5), rep("Value", 3))
  }

  # years to show in time series example plots
  to_show <- (max(burn)):(max(burn)+years_to_show)

  par(mfrow = c(10, 1), mar = c(0,0,0,0), oma = oma, cex =
    0.7, las = 1, xpd = FALSE)
  par(tck = -0.04, mgp = c(2, .45, 0))

  # environmental signal
  plot(x$env_ts[to_show], ylim = c(min(x$env_ts[to_show]),
      max(x$env_ts[to_show]) * 1.15), type = "l", xaxs = "i", axes =
    FALSE)
  box(col = "grey50")
  abline(h = 0, lty = 2, lwd = 1.5, col = "grey50")
  annotate("Environmental signal", adj = adj)
  my.axis(2, shade_years = shade_years, ylab = ylabs[1], yticks = yticks[[1]])

  # Ricker a values:
  plot_panel_lines(x$A_params[to_show, ], col = pal, ystretch = 1.3)
  my.axis(2, shade_years = shade_years, ylab = ylabs[2], yticks = yticks[[2]])
  annotate("Productivity parameter", adj = adj)

  # returns:
  plot_panel_lines(x$A[to_show, ], col = pal)
  annotate("Returns", adj = adj)
  my.axis(2, shade_years = shade_years, ylab = ylabs[3], yticks = yticks[[3]])

  # catch
  plot_panel_lines(x$F[to_show, ], col = pal)
  my.axis(2, shade_years = shade_years, ylab = ylabs[4], yticks = yticks[[4]])
  annotate("Fisheries catch", adj = adj)

  # escapement
  plot_panel_lines(x$E[to_show, ], col = pal)
  my.axis(2, shade_years = shade_years, ylab = ylabs[5], yticks = yticks[[5]])
  annotate("Escapement", adj = adj)

  # strays leaving
  plot_panel_lines(x$Strays_leaving[to_show, ], col = pal)
  my.axis(2, shade_years = shade_years, ylab = ylabs[6], yticks = yticks[[6]])
  annotate("Strays leaving", adj = adj)

  # strays joining
  plot_panel_lines(x$Strays_joining[to_show, ], col = pal)
  my.axis(2, shade_years = shade_years, ylab = ylabs[7], yticks = yticks[[7]])
  annotate("Strays joining", adj = adj)

  # SR residuals
  plot_panel_lines(x$Eps[to_show, ], ymin = "min", col = pal, ystretch = 1.4)
  abline(h = 0, lty = 2, lwd = 1.5, col = "grey50")
  my.axis(2, shade_years = shade_years, ylab = ylabs[8], yticks = yticks[[8]])
  annotate("Spawner-return residuals", adj = adj)

  plot_panel_lines(x$Est_a[to_show, ], ymin = "min", col = pal)
  annotate("Estimated a", adj = adj)
  my.axis(2, shade_years = shade_years, ylab = ylabs[9], yticks = yticks[[9]])

  plot_panel_lines(x$Est_b[to_show, ], ymin = "min", col = pal, ystretch = 1.4)
  annotate("Estimated b", adj = adj)
  my.axis(2, shade_years = shade_years, ylab = ylabs[10], yticks = yticks[[10]])

  axis(1, col = "grey50")
  mtext("Generation", side = 1, outer = TRUE, line = 2.0, cex = 0.8)
}
