"ascgen" <- function(xy = NULL, cellsize = NULL,
                     nrcol = 10, count = TRUE)
  {
      ## Verifications
      if (ncol(xy)!=2)
          stop("xy should have two columns")

      ## Remove the possible missing values
      xy <- xy[!is.na(xy[,1]),]
      xy <- xy[!is.na(xy[,2]),]


      ## Identifies the axis on which the points cover the maximum range
      xl<-c(min(xy[,1]), max(xy[,1]))
      yl<-c(min(xy[,2]), max(xy[,2]))
      rx<-xl[2]-xl[1]
      ry<-yl[2]-yl[1]
      u<-rx
      ref<-"x"
      if (ry>rx) {
          u<-ry
          ref<-"y"
      }

      ## xll and yll attributes
      xll<-xl[1]
      yll<-yl[1]

      if (!is.null(cellsize)) {
          cx<-ceiling(rx/cellsize)+1
          cy<-ceiling(ry/cellsize)+1
          asc<-matrix(0, nrow=cx, ncol=cy)
          attr(asc, "xll")<-xll
          attr(asc, "yll")<-yll
          attr(asc, "cellsize")<-cellsize
          attr(asc, "type")<-"numeric"
          class(asc)<-"asc"
      } else {
          asc<-matrix(0, nrow=nrcol, ncol=nrcol)
          cellsize<-u/(nrcol-1)
          attr(asc, "xll")<-xll
          attr(asc, "yll")<-yll
          attr(asc, "cellsize")<-cellsize
          attr(asc, "type")<-"numeric"
          class(asc)<-"asc"
      }

      ## If count TRUE, the number of points is added for each pixel
      if (count) {
          kasc<-as.kasc(list(a=asc))
          asc<-count.points(xy, kasc)
      }

      ## Output
      return(asc)
  }

