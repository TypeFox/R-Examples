################################################################################
# File: depth.graph.r
# Created by:       Oleksii Pokotylo
# First published:  01.10.2014
# Last revised:     01.10.2014
# 
# Builds the data depth surfaces for 2-dimensional data  
################################################################################

depth.graph <- function (data, depth_f = c("halfspace", "Mahalanobis", "projection", "simplicial", "simplicialVolume", "spatial", "zonoid", "none"), apoint = NULL
                           , main = depth_f
                           , xlim = c(min(data[,1]), max(data[,1])), ylim = c(min(data[,2]), max(data[,2])), zlim = c(0,max(z))
                           , xnum = 250, ynum = 250
                           , theta=15, phi=60, bold = F, ...){
  
  x1 <- seq(xlim[1], xlim[2], length = xnum)
  x2 <- seq(ylim[1], ylim[2], length = ynum)
  x1.step <- (x1[2]-x1[1])
  x2.step <- (x2[2]-x2[1])
  all.points <- as.matrix(expand.grid(x1, x2))
  all.depths <- rep(0, nrow(all.points))
  #library(depth)
  df = depth_f
  if (!is.function(depth_f)){
    depth_f = match.arg (depth_f)            
    df = switch(depth_f,
                "none" = function(x, X,...) (0),
                "zonoid" = depth.zonoid,
                "halfspace" = depth.halfspace,
                "simplicialVolume" = depth.simplicialVolume,
                "simplicial" = depth.simplicial,
                "Mahalanobis" = function(x, X,...) (.Mahalanobis_depth(x, colMeans(X), solve(cov(X)))),
                "projection" = depth.projection,
                "spatial" = depth.spatial
    )
    if (depth_f == "none") zlim = c(0,1)
  }
  
  all.depths = df(all.points, data[,1:2], ...)
  z <- matrix(all.depths, ncol=ynum, nrow=xnum, byrow=FALSE)
  
  z.red <- as.integer((data[,1]-x1[1])/x1.step+1) + as.integer((data[,2]-x2[1])/x2.step+1)*(xnum-1)
  if (bold)
    z.red <- c(z.red, 
               as.integer((data[,1]-x1[1])/x1.step+2) + as.integer((data[,2]-x2[1])/x2.step+1)*(xnum-1),
               as.integer((data[,1]-x1[1])/x1.step+1) + as.integer((data[,2]-x2[1])/x2.step+2)*(xnum-1),
               as.integer((data[,1]-x1[1])/x1.step+0) + as.integer((data[,2]-x2[1])/x2.step+1)*(xnum-1),
               as.integer((data[,1]-x1[1])/x1.step+1) + as.integer((data[,2]-x2[1])/x2.step+0)*(xnum-1)
               )
  z.black <- ifelse (is.null(apoint) || !is.numeric(apoint) || length(apoint) != 2, NA,
                     as.integer((apoint[1]-x1[1])/x1.step+1) + as.integer((apoint[1]-x2[1])/x2.step+1)*(xnum-1))
  
  zfacet <- z[-1, -1] + z[-1, -ynum] + z[-xnum, -1] + z[-xnum, -ynum]
  z.indices.zero <- which(zfacet == 0)
  cols <- rep("gray", (xnum-1)*(ynum-1))
  cols <- replace(cols, z.indices.zero, ifelse (depth_f == "none", NA,"lightblue"))
  cols <- replace(cols, z.red, "red")
  cols <- replace(cols, z.black, "black")
  
  par(bg = "white")
  persp(x1, x2, z, xlim=xlim, ylim=ylim, zlim=zlim, r = 10, theta=theta, phi=phi, 
        col=cols, main = main,
        ltheta=55, shade=0.55, ticktype="detailed", 
        xlab="x", ylab="y", zlab="D(x|X)", border=NA, box=FALSE, ...)
}
