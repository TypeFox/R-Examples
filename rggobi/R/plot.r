
# Set plot range
# Set the range of a plot.
#
# Only applies to scatterplots.
#
# @keyword internal
# @arguments GGobiDisplay object
# @arguments plot number
# @arguments x range
# @arguments y range
#X g <- ggobi(mtcars)
#X d <- displays(g)[[1]]
#X ggobi_display_set_plot_range(d, x=c(0, 40), y=c(0, 100))
# ggobi_display_set_plot_range <- function (gd, plot=1, x, y) {
#   x <- rep(x, length=2)
#   y <- rep(y, length=2)
# 
#   .GGobiCall("setPlotRange", as.numeric(x)[1], as.numeric(y)[1], as.numeric(x)[2], as.numeric(y)[2], gd, as.integer(plot))
#   invisible()
# }

# Get plot range
# Get plot range as list of X and Y ranges
# 
# @keyword internal
# @arguments GGobiDisplay object
# @arguments plot number
#X g <- ggobi(mtcars)
#X d <- displays(g)[[1]]
#X ggobi_display_get_plot_range(d)
# ggobi_display_get_plot_range <- function (gd, plot=1) {
#   range <- .GGobiCall("getPlotRange", gd, as.integer(plot))
#   names(range) <- c("x", "y")
#   range
# }
