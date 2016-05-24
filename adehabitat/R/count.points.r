"count.points" <- function(xy, w)
  {
      ## Verifications
      if (inherits(w, "asc"))
          w<-as.kasc(list(toto=w))
      if (inherits(w, "kasc"))
          w<-storemapattr(w)
      if (!inherits(w, "mapattr"))
          stop("non convenient format for w")

      ## Prepares a vector containing the boundaries of the pixels
      xyc<-getXYcoords(w)
      xc<-xyc$x-attr(w, "cellsize")/2
      yc<-xyc$y-attr(w, "cellsize")/2
      xc<-c(xc, max(xc)+attr(w, "cellsize"))
      yc<-c(yc, max(yc)+attr(w, "cellsize"))

      ## discretize the points according to these classes
      x<-xy[,1]
      y<-xy[,2]
      x<-cut(x, xc)
      y<-cut(y, yc)

      ## Transform into an object of class "asc"
      output<-as.matrix(table(x, y))
      attr(output, "xll")<-attr(w, "xll")
      attr(output, "yll")<-attr(w, "yll")
      attr(output, "cellsize")<-attr(w, "cellsize")
      attr(output, "type")<-"numeric"
      class(output)<-"asc"
      return(output)
  }

