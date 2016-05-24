#' Plot Profile Intersection
#' 
#' Plots resistivity against depth on and next to the intersection line 
#' between two profiles.
#' 
#' @param .Object1 either a single Profile or a ProfileSet.
#' @param .Object2 either a second single Profile or NULL if .Object1 is of type ProfileSet.
#' @param main title to be plotted.
#' @param xlab label of the x-axes, e.g. length [m].
#' @param ylab label of the y-axes, e.g. height above sea level [m].
#' @param col character vector of colors.
#' @param pch numeric vector of plotting symbols.
#' @param type plot type (default "p" for points).
#' "b" for both points and lines, "c" for empty points joined by lines, 
#' "o" for overplotted points and lines, 
#' "s" and "S" for stair steps and "h" for histogram-like vertical lines. 
#' Finally, "n" does not produce any points or lines.
#' @param legendLoc legendLocation (default "bottomleft").
#' @param trafo transformation to be done on data (default: log).
#' @param backtrafo back transformation to plot correct labels (default: exp).
#' @export
#' @seealso \code{\link{ProfileSet-class}}
#' @examples 
#' # data(sinkhole)
#' 
#' # plotIntersect(sinkhole)
#' # plotIntersect(sinkhole@profiles[[1]], sinkhole@profiles[[2]])
setGeneric("plotIntersect", function(.Object1, .Object2=NULL, 
                                     xlab="Height above sea level [m]", 
                                     ylab=expression(paste("Resistivity [", Omega, "m]")), 
                                     main="", trafo=log, backtrafo=exp,
                                     col=colors, pch=c(20,20), type="p",
                                     legendLoc="bottomleft") {
  standardGeneric("plotIntersect")  
})

#' @rdname plotIntersect
#' @export
setMethod("plotIntersect", signature(.Object1="ProfileSet"),
          function(.Object1, xlab, ylab, main, trafo, backtrafo, col, pch, type, 
                   legendLoc) {
            for(i in 1:(length(.Object1@profiles)-1)) 
              for(j in (i+1):length(.Object1@profiles))
                plotIntersect(.Object1@profiles[[i]], .Object1@profiles[[j]],
                              xlab, ylab, main, trafo, backtrafo, col, pch, type, legendLoc)
          })

#' @rdname plotIntersect
#' @export
setMethod("plotIntersect", signature(.Object1="Profile", .Object2="Profile"),
          function(.Object1, .Object2, xlab, ylab, main, trafo, backtrafo, 
                   col, pch, type, legendLoc) {
            # slopes m
            m1 <- .Object1@gpsCoordinates@lmRelative$coefficients[2]
            m2 <- .Object2@gpsCoordinates@lmRelative$coefficients[2]
            
            # intercepts n
            n1 <- .Object1@gpsCoordinates@lmRelative$coefficients[1]
            n2 <- .Object2@gpsCoordinates@lmRelative$coefficients[1]
            
            # calculate intersection point
            # m1 * x.intersect + n1 = m2 * x.intersect + n2
            x.intersect <- (n2 - n1)/(m1 - m2)
            y.intersect <- m1 * x.intersect + n1
            
            # starting points of Profile 1 and 2
            x.start1 <- min(.Object1@gpsCoordinates@relative$lon)
            y.start1 <- m1 * x.start1 + n1            
            x.start2 <- min(.Object2@gpsCoordinates@relative$lon)
            y.start2 <- m2 * x.start2 + n2
            
            # calculate length (hypotenuse) from starting to intersection point
            x.diff1 <- x.intersect - x.start1
            y.diff1 <- y.intersect - y.start1
            x.diff2 <- x.intersect - x.start2
            y.diff2 <- y.intersect - y.start2
            
            length1 <- sqrt(x.diff1^2 + y.diff1^2)
            length2 <- sqrt(x.diff2^2 + y.diff2^2)
            
            # identify point indices on intersection line and next to it            
            indices1 <- c(which(round(.Object1@xyzData@heightAdaption$dist) == round(length1)),
                          which(round(.Object1@xyzData@heightAdaption$dist) == round(length1 + 1)),
                          which(round(.Object1@xyzData@heightAdaption$dist) == round(length1 - 1)))
            
            indices2 <- c(which(round(.Object2@xyzData@heightAdaption$dist) == round(length2, 0)),
                          which(round(.Object2@xyzData@heightAdaption$dist) == round(length2 + 1, 0)),
                          which(round(.Object2@xyzData@heightAdaption$dist) == round(length2 - 1, 0)))
            
            # check whether there is an intersection
            if(length(indices1) == 0 |length(indices2) == 0) {
              print(paste("No intersection between ", .Object1@title, " and ", .Object2@title, ".", sep=""))
              return()
            }
            
            # identify xyz values for these indices
            res1 <- data.frame(
              "dist" = .Object1@xyzData@heightAdaption$dist[indices1],
              "depth" = .Object1@xyzData@heightAdaption$depth[indices1],
              "val" = .Object1@xyzData@heightAdaption$val[indices1])
            
            res2 <- data.frame(
              "dist" = .Object2@xyzData@heightAdaption$dist[indices2],
              "depth" = .Object2@xyzData@heightAdaption$depth[indices2],
              "val" = .Object2@xyzData@heightAdaption$val[indices2])
            
            #boxplot(trafo(res1$val)~round(res1$depth))
            lab.breaks <- round(backtrafo(seq(trafo(min(res1$val, res2$val)),
                                      trafo(max(res1$val, res2$val)), 
                                      length.out=6)))
            at.breaks <- seq(trafo(min(res1$val, res2$val)),
                                      trafo(max(res1$val, res2$val)), 
                                      length.out=6)
            
            plot(res1$depth, trafo(res1$val),
                 xlim=c(min(res1$depth, res2$depth), max(res1$depth, res2$depth)),
                 ylim=c(trafo(min(res1$val, res2$val)), trafo(max(res1$val, res2$val))),
                 xlab=xlab, ylab=ylab, main=main, col=col[1], pch=pch[1], type=type,
                 yaxt="n")
            points(res2$depth, trafo(res2$val), col=col[2], pch=pch[2], type=type)   
            legend(legendLoc, col=col, pch=pch, 
                   legend=c(.Object1@title, .Object2@title))
            axis(side=2, at=at.breaks, labels=lab.breaks)
          })