#' Plot Raw Data Points
#' 
#' Plots raw data points of a single profile (e.g. to show measurement gaps). 
#' @param Profile profile.
#' @param xlab label for x-axes.
#' @param ylab label for y-axes.
#' @param main title to be plotted.
#' This can either be a single character or an integer code for one of a set of graphics symbols.
#' @param ... plot arguments (like pch, cex, col,...).
#' @export
#' @seealso \code{\link{plotRawHeight}}, \code{\link{RawData-class}}
#' @examples
#' data(sinkhole) 
#' plotRaw(sinkhole@profiles[[1]])
#' plotRaw(sinkhole@profiles[[2]])
#' plotRaw(sinkhole@profiles[[3]])
plotRaw <- function(Profile, xlab="Length [m]", ylab="Depth [m]",
                    main=paste(Profile@title, "without topography"), ...) {  
  plot(Profile@rawData@seaLevel$dist, -1*(Profile@rawData@seaLevel$depth), 
       xlab=xlab, ylab=ylab, main=main, asp=1, ...) 
}

#' Plot Raw Data Points considering Topography
#' 
#' Plots raw data points of a single profile (e.g. to show measurement gaps).
#' The topography is considered, i.e., heights are added to the measurement depth. 
#' @param Profile profile.
#' @param height topo data frame of distances and height.
#' Use "Profile@xyzData@height" instead.
#' @param spline if TRUE spline interpolation is conducted.
#' @param xlab label for x-axes.
#' @param ylab label for y-axes.
#' @param main title to be plotted.
#' @param ... plot arguments (like pch, cex, col,...).
#' @export
#' @seealso \code{\link{plotRaw}}, \code{\link{RawData-class}}
#' @examples 
#' data(sinkhole)
#' plotRawHeight(sinkhole@profiles[[2]])
#' plotRawHeight(sinkhole@profiles[[2]], sinkhole@profiles[[2]]@xyzData@height)
plotRawHeight <- function(Profile, height=Profile@rawData@height, 
                          spline=TRUE,
                          xlab="Length [m]", ylab="Depth [m]",
                          main=paste(Profile@title, "without topography"), ...) { 
  if(spline) {
  height <- as.data.frame(spline(height, method="natural",
                   xmin=min(Profile@rawData@seaLevel$dist), 
                   xmax=max(Profile@rawData@seaLevel$dist),
                   n=(max(Profile@rawData@seaLevel$dist)-min(Profile@rawData@seaLevel$dist)+1)))
  }

  heightAdaption <- data.frame(
    "dist"=Profile@rawData@seaLevel$dist, 
    "depth"=-1*Profile@rawData@seaLevel$depth, 
    "val"=Profile@rawData@seaLevel$val)

  for(i in 1:nrow(height)) {
    indices <- which(round(height[i,1]) == round(Profile@rawData@seaLevel$dist))
    if (length(indices) > 0)
      for(j in 1:length(indices)) {
        heightAdaption$depth[indices[j]] <- heightAdaption$depth[indices[j]] + height[i,2] 
      }
  }         
  plot(heightAdaption$dist, heightAdaption$depth, 
       xlab=xlab, ylab=ylab, main=main, asp=1, ...) 
}