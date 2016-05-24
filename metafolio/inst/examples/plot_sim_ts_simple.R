# plot time series from the simulation
# just the annual rate of growth and returns in this one

plot_sim_ts_simple <- function(x, pal, burn = 1:60, text = "Returns") {
# x is output from a simulation run; it's a list
# pal is the colour palette

to_show <- (max(burn)):(max(burn)+40) # years to show in time series example plots

annotate <- function(label, xfrac = 0.008, yfrac = 0.18, pos = 4, cex = 0.9, ...) {
  u <- par("usr")
  x <- u[1] + xfrac * (u[2] - u[1])
  y <- u[4] - yfrac * (u[4] - u[3])
  text(x, y, label, pos = pos, cex = cex, ...)
}

#my.axis <- function(side) axis(side, col = "grey50", at = c(800, 1200))

matplot(x$A[to_show, 1:3], type = "l", col = pal, lty = 1, ylab = "Returns", xlab = "Time", xaxt = "n", axes = FALSE, xaxs = "i", ylim = c(0, 6500), yaxs = "i")
annotate(text)
box(col = "grey50")
axis(2, at = c(0, 5000), col = "grey50", tck = -0.05, mgp = c(2, 0.5, 0))
#my.axis(1)

}

