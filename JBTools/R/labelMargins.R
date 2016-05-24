labelMargins <- function(
  ##title<< Label the margins of a plot
  labels     ##<< character vector: labels to use
  , side = 1 ##<< integer: side of the plot to label
  , ...      ##<< further arguments passed to the plotting routines
  )
  ##description<< Writes equidistant labels to the margins of a plot.
{
  at <- 1/length(labels) * 0:(length(labels) - 1)+ 1/(2*length(labels))
  if (side == 2)
    at = rev(at)
  mtext(labels, at = at, side = side, ...)
  ##value<< nothing is returned
}
