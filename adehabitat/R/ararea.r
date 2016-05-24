"ararea" <- function(x)
  {
      ## Verifications
      if (!inherits(x, "area"))
          stop("x should be of class \"area\"")

      ## Computes the area of each polygon
      uu <- split(x[,2:3], x[,1])
      foo <- function(y) {
          class(y) <- "data.frame"
          if (!all(unlist(y[1,])==unlist(y[nrow(y),])))
              y <- rbind(y,y[1,])
          pol <- Polygon(as.matrix(y))
          spdf <- SpatialPolygons(list(Polygons(list(pol), 1)))
          lar <- unlist(lapply(polygons(spdf)@polygons,
                               function(x) unlist(lapply(x@Polygons, function(y)
                                                         .arcp(y@coords)))))
          lhol <- unlist(lapply(polygons(spdf)@polygons,
                                function(x) unlist(lapply(x@Polygons, function(y)
                                                          y@hole))))
          sum(lar[!lhol])-sum(lar[lhol])
      }

      ## Output
      res <- unlist(lapply(uu, foo))
      names(res) <- names(uu)
      return(res)
  }

