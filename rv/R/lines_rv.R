
# ========================================================================
# lines.rv  -  plot some random lines
# ========================================================================
# btw, "rvlines" does not make sense - that'd be something like a weird histogram

lines.rv <- function(x, y, type="l", ...) {
  if (is.rvobj(x) || is.rvobj(y)) {
    points.rv(x, y, type="l", ...)
  } else {
    lines(x, y, type=type, ...)
  }
}

