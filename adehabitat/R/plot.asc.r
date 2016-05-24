"plot.asc" <- function(x, ...)
  {
      ## Verifications
      if (!inherits(x, "asc"))
          stop("should be an object of class \"asc\"")
      if (attr(x, "type")=="factor")
          stop("not implemented for factors")

      ## The plot
      xy<-getXYcoords(x)
      filled.contour(xy$x, xy$y, x, asp=1, ...)
  }

