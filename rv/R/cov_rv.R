# ========================================================================
# cov  -  short description
# ========================================================================


cov.rv <- function(x, y=NULL, ...)  ## EXPORT cov.rv
{
  if (!is.matrix(x)) {
    if (is.null(y)) {
      stop("supply both x and y or a matrix-like x")
    }
    x <- as.vector(x)
  }
  if (is.rvobj(x) || is.rvobj(y)) {
    rvmapply(stats:::cov, x=x, y=y, ...)
  } else {
    stats:::cov(x=x, y=y, ...)
  }
}



