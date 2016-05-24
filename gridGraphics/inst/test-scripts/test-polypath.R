
library(gridGraphics)

plotPath <- function(x, y, col = "grey", rule = "winding") {
    plot.new()
    plot.window(range(x, na.rm = TRUE), range(y, na.rm = TRUE))
    polypath(x, y, col = col, rule = rule)
    if (!is.na(col))
        mtext(paste("Rule:", rule), side = 1, line = 0)
}

plotRules <- function(x, y, title) {
    plotPath(x, y)
    plotPath(x, y, rule = "evenodd")
    mtext(title, side = 3, line = 0)
    plotPath(x, y, col = NA)
}

polypath1 <- function() {
    par(mfrow = c(5, 3), mar = c(2, 1, 1, 1))
    plotRules(c(.1, .1, .9, .9, NA, .2, .2, .8, .8),
              c(.1, .9, .9, .1, NA, .2, .8, .8, .2),
              "Nested rectangles, both clockwise")
    plotRules(c(.1, .1, .9, .9, NA, .2, .8, .8, .2),
              c(.1, .9, .9, .1, NA, .2, .2, .8, .8),
              "Nested rectangles, outer clockwise, inner anti-clockwise")
    plotRules(c(.1, .1, .4, .4, NA, .6, .9, .9, .6),
              c(.1, .4, .4, .1, NA, .6, .6, .9, .9),
              "Disjoint rectangles")
    plotRules(c(.1, .1, .6, .6, NA, .4, .4, .9, .9),
              c(.1, .6, .6, .1, NA, .4, .9, .9, .4),
              "Overlapping rectangles, both clockwise")
    plotRules(c(.1, .1, .6, .6, NA, .4, .9, .9, .4),
              c(.1, .6, .6, .1, NA, .4, .4, .9, .9),
              "Overlapping rectangles, one clockwise, other anti-clockwise")
}

plotdiff(expression(polypath1()), "polypath-1")

plotdiffResult()

