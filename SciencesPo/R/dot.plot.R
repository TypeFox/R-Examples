#' @encoding UTF-8
#' @title Dot Plot
#'
#' @description Makes a dot plot.
#' @param x The data vector
#' @param pch The plotting "character" or symbol, default is dots.
#' @param bins The bins width.
#' @param spacing A value for vertically spacing between dots.
#' @param xlab The axis label.
#' @param \dots Other parameters passed on to `plot`.
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#' @export
#' @examples
#' with(iris, dot.plot(Sepal.Length,
#' xlab="Length", col = as.numeric(Species)))
#'
`dot.plot` = function(x, pch = 16, bins = 50, spacing = 1, xlab, ...) {
  if(missing(xlab))
    xlab = as.character(substitute(x))

  # determine dot positions
  inc = diff(pretty(x, n = bins)[1:2])
  freq = table(inc * round(x / inc, 0))
  xx = rep(as.numeric(names(freq)), freq)
  yy = unlist(lapply(freq, seq_len))

  # make the order of the dots the same as the order of the data
  idx = seq_along(x)
  idx[order(x)] = idx
  xx = xx[idx]
  yy = yy[idx]

  # make a blank plot
  graphics::plot(xx, yy, type = "n", axes = FALSE, xlab = xlab, ylab = "")

  # draw scale
  graphics::axis(1)
  ylow = graphics::par("usr")[3]
  graphics::abline(h = ylow) # extend to full width

  # draw points and support resizing
  grDevices::recordGraphics({
    yinc = 0.5 * spacing * graphics::par("cxy")[2]
    graphics::points(xx, ylow + yinc * (yy - .5), pch = pch, ...)
  },
  list(),
  environment(NULL))

  invisible()
}
NULL
