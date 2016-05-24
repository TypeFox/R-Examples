
# Function to plot the x (time) axis for overlapPlot and densityPlot

# Not exported

plotTimeAxis <- function(xscale, ...) {
  # Deal with ... argument:
  dots <- list(...)
  if(length(dots) == 1 && class(dots[[1]]) == "list")
    dots <- dots[[1]]
  selPlot <- names(dots) %in%
     c("cex.axis", "col.axis", "family", "font.axis", "las", "tck", "tcl", "xaxt",
        "tick", "lwd.ticks", "col.ticks")
  plotArgs <- dots[selPlot]
  plotArgs$side <- 1

  if(is.na(xscale)) {
      plotArgs$at <- c(-pi, -pi/2, 0, pi/2, pi, 3*pi/2, 2*pi)
      plotArgs$labels <- c(expression(-pi), expression(-pi/2), "0",
          expression(pi/2), expression(pi),
          expression(3*pi/2), expression(2*pi))
  } else if(xscale == 24) {
    plotArgs$at <- c(-12, -6, 0,6,12,18,24)
    plotArgs$labels <- c("12:00", "18:00", "0:00", "6:00", "12:00", "18:00", "24:00")
  } else if(xscale == 1) {
    plotArgs$at <- c(-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)
    plotArgs$labels <- TRUE
  } 
  do.call(axis, plotArgs, quote=TRUE)
}