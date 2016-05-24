AddSegments <- function (x, y, half.width.factor = 0.45, ...) {
  ## Add horizontal line segments at the specified x and y points.
  ## The main use for this function is to add reference lines to
  ## side-by-side boxplots.
  ##
  ## Args:
  ##   x: A vector of equally spaced values on which the segments
  ##     should center.
  ##   y:  The value of each line segment on the vertical axis.
  ##   half.width.factor: What fraction of the distance between
  ##     consecutive x values should be taken up by the line segments.
  ##     .5 == all of it.
  ##   ...: Extra arguments to be passed to 'segments'
  ##
  dx = diff(x)
  if (!(all(dx) == dx[1])) {
    stop("need equally spaced x's in add.segments")
  }
  dx <- dx[1]
  half.width = half.width.factor * dx
  x0 <- x - half.width
  x1 <- x + half.width
  segments(x0, y, x1, y, ...)
}
