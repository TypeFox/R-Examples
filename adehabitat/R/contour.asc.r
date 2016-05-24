"contour.asc" <- function(x, ...)
  {
      ## Verifications
      if (!inherits(x, "asc")) stop("not an \"asc\" object")
      if (attr(x, "type")=="factor")
          stop("function contour cannot be used with factors")

      ## Use of the function contour
      z<-x
      xy<-getXYcoords(z)
      x<-xy$x
      y<-xy$y
      contour(x=x, y=y, z,  ...)
  }

