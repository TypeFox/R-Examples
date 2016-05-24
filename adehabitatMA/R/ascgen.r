"ascgen" <- function(xy, cellsize = NULL,
                     nrcol = NULL, count = TRUE)
  {
      ## Verifications
      if (!inherits(xy, "SpatialPoints"))
          stop("xy should inherit the class SpatialPoints")
      if (is.null(cellsize)&is.null(nrcol))
          stop("One of the parameters cellsize or nrcol should be specified")
      if (!is.null(cellsize)&!is.null(nrcol))
          warning("cellsize and nrcol specified.\nOnly cellsize is taken into account")
      if (ncol(coordinates(xy))>2)
          stop("xy should be defined in two dimensions")

      pxy <- proj4string(xy)

      ## Remove the possible missing values
      xy <- coordinates(xy)

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

      if (is.null(cellsize)) {
          cellsize <- u/nrcol
      }

      cx<-ceiling(rx/cellsize)+1
      cy<-ceiling(ry/cellsize)+1

      xx <- seq(xl[1], xl[2]+cellsize, by=cellsize)
      yy <- seq(yl[1], yl[2]+cellsize, by=cellsize)
      grid <- expand.grid(yy,xx)[,2:1]
      asc<-data.frame(x=rep(0, nrow(grid)))
      coordinates(asc) <- grid
      gridded(asc) <- TRUE

      ## If count TRUE, the number of points is added for each pixel
      if (count) {
          uu <- over(SpatialPoints(xy), geometry(asc))
          uu <- table(uu)
          asc <- asc[[1]]
          asc[as.numeric(names(uu))] <- uu
          asc <- data.frame(x=asc)
          coordinates(asc) <- grid
          gridded(asc) <- TRUE
      }
      if (!is.na(pxy))
          proj4string(asc) <- CRS(pxy)

      ## Output
      return(asc)
  }

