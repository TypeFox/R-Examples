plot.gambin <-
function(x, barcol = "grey", barwidth = 1, cex.dots = 1, dotpch = 16, dotcol = par("fg"), line = FALSE, lwd = 1, linecol = par("fg"), ...)
{
  ylim <- max(predict(x), x$Data$species)
  midpoints <- barplot(x$Data$species, names.arg = x$Data$octave, ylim = c(0, 0.5 + ylim), col = barcol, width = barwidth,  xlab = "Octaves", ylab = "Number of species", ...)
  points(midpoints, predict(x), pch = dotpch, cex = cex.dots, col = dotcol)
  if(line)
      lines(midpoints, predict(x), lwd = lwd, col = linecol)
}
