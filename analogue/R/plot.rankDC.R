plot.rankDC <- function(x, sort = TRUE, decreasing = FALSE,
                        xlab = "Rank correlation", color = "blue",
                        pch = 21, bg = "blue", lcolor = "grey",
                        lty = "solid", ...) {
  if(sort)
    x <- sort(x, decreasing = decreasing)
  dotchart(x, xlab = xlab,
           color = color, pch = pch, bg = bg, lcolor = lcolor,
           lty = lty, ...)
}
