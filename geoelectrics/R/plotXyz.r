#' Plot XYZ Points
#' 
#' Plots the interpolated points of the xyz data. 
#' @param Profile profile.
#' @param xlab label for x-axes.
#' @param ylab label for y-axes.
#' @param main title to be plotted.
#' @param asp the y/x aspect ratio (default: 1).
#' @param ... plot parameters
#' @export
#' @seealso \code{\link{Profile-class}}, \code{\link{plotXyzHeight}}, 
#' \code{\link{plot3dXyz}}, 
#' \code{\link{levelplotXyz}}, \code{\link{levelplotXyzHeight}}
#' @examples 
#' data(sinkhole)
#' 
#' plotXyz(sinkhole@profiles[[1]])
plotXyz <- function(Profile, xlab="Length [m]", ylab="Depth [m]",
                    main=paste(Profile@title, "without topography"),
                    asp=1, ...) {
  plot(Profile@xyzData@seaLevel$dist, Profile@xyzData@seaLevel$depth, 
       xlab=xlab, ylab=ylab, main=main, asp=asp, ...)
}

#' Plot XYZ Points considering Topography
#' 
#' Plots the interpolated points of the xyz data
#' after height adjustment. 
#' @param Profile profile.
#' @param xlab label for x-axes.
#' @param ylab label for y-axes.
#' @param main title to be plotted.
#' @param asp the y/x aspect ratio (default: 1).
#' @param ... plot parameters.
#' @export
#' @seealso \code{\link{Profile-class}}, \code{\link{plotXyz}}, 
#' \code{\link{plot3dXyz}}, 
#' \code{\link{levelplotXyz}}, \code{\link{levelplotXyzHeight}}
#' @examples 
#' data(sinkhole)
#' 
#' plotXyzHeight(sinkhole@profiles[[1]])
plotXyzHeight <- function(Profile, xlab="Length [m]", ylab="Depth [m]",
                          main=paste(Profile@title, "with topography"), asp=1, ...) {
  plot(Profile@xyzData@heightAdaption$dist, Profile@xyzData@heightAdaption$depth, 
       xlab=xlab, ylab=ylab, main=main, asp=asp, ...)
}