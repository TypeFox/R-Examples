#' Plots Profiles 3D
#' 
#' Plots the interpolated resistance values of the 
#' xyz data for all profiles.
#' 
#' @param .Object either an object of a single Profile or a ProfileSet.
#' @param title title to be plotted.
#' @param sub subtitle to be plotted.
#' @param xlab label of the x-axes, e.g. length [m].
#' @param ylab label of the y-axes, e.g. height above sea level [m].
#' @param zlab label of the z-axes, e.g. length [m].
#' @param minData mimimum value to adjust color bar.
#' @param maxData maximum value to adjust color bar.
#' @param col vector of colors.
#' @param trafo transformation to be done on data (default: log).
#' @param psize size of value points (default: 10).
#' @export
#' @seealso \code{\link{Profile-class}}, \code{\link{ProfileSet-class}},
#' \code{\link{plotXyz}}, \code{\link{levelplotXyz}}
#' @examples
#' # data(sinkhole)
#' 
#' # plot3dXyz(sinkhole@profiles[[1]])
#' # plot3dXyz(sinkhole)
setGeneric("plot3dXyz", function(.Object, title="", sub="",
                                 xlab="", ylab="", zlab="",
                                 minData=0, maxData=9999999, 
                                 col=colors, trafo=log, psize=pointsize){
  standardGeneric("plot3dXyz")
})

#' @rdname plot3dXyz
#' @aliases plot3d
#' @export
setMethod("plot3dXyz", signature(.Object="ProfileSet"),
          function(.Object, title=.Object@title, sub="",
                   xlab="", ylab="", zlab="",
                   minData=.Object@minData,
                   maxData=.Object@maxData, col, trafo, psize) {
            lapply(.Object@profiles, plot3dXyz,
                   minData=minData, maxData=maxData, col=col, trafo=trafo)
            title3d(title, sub, xlab, ylab, zlab)
          })

#' @rdname plot3dXyz
#' @aliases plot3d
#' @export
setMethod("plot3dXyz", signature(.Object="Profile"),
          function(.Object, title="", sub="",
                   xlab="", ylab="", zlab="",
                   minData=.Object@xyzData@minData, 
                   maxData=.Object@xyzData@maxData, 
                   col, trafo, psize) {
            title3d(title, sub, xlab, ylab, zlab)
            values <- trafo(.Object@xyzData@heightAdaption$val)
            colorAssignment <- myColorRamp(col, values, trafo(minData), trafo(maxData))
            
            l <- .Object@xyzData@heightAdaption$dist # hypotenuse
            m <- .Object@gpsCoordinates@lmRelative$coefficients[2] # y = mx + n
            n <- .Object@gpsCoordinates@lmRelative$coefficients[1]
            alpha <- atan(m)
            
            # calculate adjacent leg
            x <- cos(alpha) * l
            
            # get starting point and adjust
            start.x <- min(.Object@gpsCoordinates@relative$lon)
            x <- x + start.x
            
            # calculate opposite leg
            y <- m * x + n
            
            # plot 3D    
            rgl.bg(color="white")
            points3d(y, .Object@xyzData@heightAdaption$depth, x, color=colorAssignment, size=pointsize)  
            rgl.bbox()  
            rgl.texts(y[1], .Object@xyzData@heightAdaption$depth[1]+20, x[1], 
                      text=paste(.Object@title), cex=1, color="black")
            axes3d(edges="bbox", yunit=25, expand=1.2)
          })