gcolors = c("red", "blue", "green", "orange", "violet")

depth.contours.ddalpha <- function(ddalpha, main = "", xlab="", ylab = "", drawplot = T, frequency=100){
  
  if(class(ddalpha)!="ddalpha")
    stop("Not a 'ddalpha' classifier")
  
  if (ddalpha$dimension != 2)
  {
    warning ("The contours may be drawn only for 2 dimensional datasets")
    return(0)
  }
  
  if (ddalpha$needtransform == 1)
    data = rbind(ddalpha$patterns[[1]]$transformer(ddalpha$patterns[[1]]$points, inv = T), ddalpha$patterns[[2]]$transformer(ddalpha$patterns[[2]]$points, inv = T))
  else
    data = rbind(ddalpha$patterns[[1]]$points, ddalpha$patterns[[2]]$points)
  classes = c(rep(ddalpha$patterns[[1]]$name, ddalpha$patterns[[1]]$cardinality), rep(ddalpha$patterns[[2]]$name, ddalpha$patterns[[2]]$cardinality))
  margins = c(min(data[,1]), max(data[,1]), min(data[,2]), max(data[,2]));  margins = margins + c(-0.1*(margins[2]-margins[1]), 0.1*(margins[2]-margins[1]), -0.1*(margins[4]-margins[3]), 0.1*(margins[4]-margins[3]))
  C = ncol(data)
  
  if (drawplot)
    plot(data, col = ifelse (classes == ddalpha$patterns[[1]]$name, "red", "blue"), main = main, xlab=xlab, ylab = ylab)
  if (is.null(ddalpha$methodDepth))
    return(0)
  
  if (ddalpha$methodDepth == "Mahalanobis" && ddalpha$mahEstimate == "moment"){
    # Mahalanobis depth
    for (i in seq(length(ddalpha$patterns))){
      mahalanobisRegions(ddalpha$patterns[[i]]$points, col = gcolors[i])
    }
  } else 
    if (ddalpha$methodDepth == "Mahalanobis" && ddalpha$mahEstimate == "MCD"){
      # Mahalanobis depth mit MCD
      for (i in seq(length(ddalpha$patterns))){
        mahalanobisMCDRegions(ddalpha$patterns[[i]]$points, ddalpha$mahParMcd, col = gcolors[i])
      }
    } else 
#       if (ddalpha$methodDepth == "zonoid"){
#         if (require(WMTregions))
#           for (i in seq(length(ddalpha$patterns))){
#             wmtRegions(ddalpha$patterns[[i]]$points, depths = c(1:9/10), trtype = "zonoid", col = gcolors[i])
#           }
#         else{
#           contourRegions(ddalpha, margins = margins, depths = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), frequency=frequency)      
#         }
#       } else 
#         if (ddalpha$methodDepth == "halfspace"){
#           library(depth)
#           for (i in seq(length(ddalpha$patterns))){
#             locationRegions(ddalpha$patterns[[i]]$points, depths = c(1:9/10), col = gcolors[i])
#           }    
#         }  else
          # no formal depth regions
        {
          return (contourRegions(ddalpha, margins = margins, depths = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), frequency=frequency))
        }   
  
  return (0)
}


depth.contours <- function(data, depth, main = "", xlab="", ylab = "", drawplot = T, frequency=100, col = "red", ...){
  
  if(!(is.matrix(data)||is.data.frame(data)))
    stop("Data is not a matrix or a data frame")
  
  if (ncol(data) != 2)
  {
    warning ("The contours may be drawn only for 2 dimensional datasets")
    return
  }
  
  if (drawplot)
    plot(data, col = col, main = main, xlab=xlab, ylab = ylab, ...)
  margins = c(min(data[,1]), max(data[,1]), min(data[,2]), max(data[,2]));  margins = margins + c(-0.1*(margins[2]-margins[1]), 0.1*(margins[2]-margins[1]), -0.1*(margins[4]-margins[3]), 0.1*(margins[4]-margins[3]))
  
  
  if (depth == "Mahalanobis"){
    # Mahalanobis depth
    mahalanobisRegions(data, col = col)
  } else 
#     if (depth == "zonoid"){
#       if (require(WMTregions))
#         wmtRegions(data, depths = c(1:9/10), trtype = "zonoid", col = col)
#       else{
#         contourRegionsData(data, depth = depth, margins = margins, depths = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), frequency=frequency, col = col, ...)      
#       }
#     } else 
#       if (depth == "halfspace"){
#         library(depth)
#         locationRegions(data, depths = c(1:9/10), col = col)  
#       }  else 
        # no formal depth regions
      {
        return (contourRegionsData(data, depth = depth, margins = margins, depths = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), frequency=frequency, col = col, ...))
      }   
  
  return(0)
}


# wmtRegions <- function(x, depths = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), trtype = "zonoid", col = "black"){
#   fname <- "Cloud.dat"
#   #trname <- "TRegion_vertices.dat"
#   trname <- "tmp_vrtheap.dat"
#   fdir = getwd()
#   for (i in 1:length(depths)){
#     ffullname = as.character(paste(fdir, "/", fname, sep = ""))
#     write(trtype, ffullname)
#     write(depths[i], ffullname, append = TRUE)
#     write(ncol(x), ffullname, append = TRUE)
#     write(nrow(x), ffullname, append = TRUE)
#     write(array(t(x), dim = c(nrow(x), ncol(x))), ffullname, ncolumns = ncol(x), append = TRUE)
#     WMTR()
#     cat("Calculated for i = ", i, "\n")
#     wmtreg <- read.table(as.character(paste(fdir, "/", trname, sep = "")), sep=" ")
#     unlink(as.character(paste(fdir, "/", trname, sep = "")))
#     wmtreg <- as.matrix(wmtreg)
#     numribs <- nrow(wmtreg)/2
#     for (i in 1:numribs){
#       lines(rbind(wmtreg[i*2 - 1,], wmtreg[i*2,]), col = col)
#     }    
#   }
# }
# 
# locationRegions <- function(x, depths = NULL, col = "black"){
#   numLearn <- nrow(x)
#   vert = isodepth(x, c(2:(numLearn - 1)), output = T, dpth = as.integer(depths*numLearn/2))
#   for (verticles in vert){
#     if(is.null(verticles) || nrow(verticles)<1) next
#     verticles = rbind(verticles, verticles[1,])
#     lines(verticles, col = col)
#   }
# }

mahalanobisRegions <- function(x, depths = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), col = "black"){
  c <- c(0,0)
  for (i in 1:nrow(x)){
    c <- c + x[i,]
  }
  mu <- c / nrow(x)
  sigma.inv <- solve(cov(x))
  for (i in 1:length(depths)){
    ellipsem(mu = mu, amat = sigma.inv, c2 = 1/depths[i] - 1, showcentre = F, col = col)
  }
}

mahalanobisMCDRegions <- function(x, alpha = 1/2, depths = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), col = "black"){
  #  library(robustbase)
  estimate <- covMcd(x, alpha = alpha)
  mu <- estimate$center
  sigma.inv <- solve(estimate$cov)
  for (i in 1:length(depths)){
    ellipsem(mu = mu, amat = sigma.inv, c2 = 1/depths[i] - 1, showcentre = F, col = col)
  }
}

contourRegions <- function(ddalpha, margins = NULL, depths = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), frequency=100){
  
  if (is.null(margins)){
    if (ddalpha$needtransform == 1)
      data = rbind(ddalpha$patterns[[1]]$transformer(ddalpha$patterns[[1]]$points, inv = T), ddalpha$patterns[[2]]$transformer(ddalpha$patterns[[2]]$points, inv = T))
    else
      data = rbind(ddalpha$patterns[[1]]$points, ddalpha$patterns[[2]]$points)
    margins = c(min(data[,1]), max(data[,1]), min(data[,2]), max(data[,2]));  margins = margins + c(-0.1*(margins[2]-margins[1]), 0.1*(margins[2]-margins[1]), -0.1*(margins[4]-margins[3]), 0.1*(margins[4]-margins[3]))
  }
  
  gx <- seq(margins[1], margins[2], length=frequency)
  gy <- seq(margins[3], margins[4], length=frequency)
  y <- as.matrix(expand.grid(gx, gy))    
  
  depthcontours <- .ddalpha.count.depths(ddalpha, y)  
  
  for (i in seq(length(ddalpha$patterns))){
    contour(gx, gy, matrix(depthcontours[,i], nrow=length(gx), ncol=length(gy)), add=TRUE, levels=depths*max(depthcontours), drawlabels=FALSE, col = gcolors[i])
  }
  
  return(depthcontours)
}

contourRegionsData <- function(data, depth, margins = NULL, depths = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), frequency=100, col = "red", ...){
  
  if (is.null(depth) || depth == "none")
    return (0);
  df = switch(depth,
              "zonoid" = depth.zonoid,
              "halfspace" = depth.halfspace,
              "simplicialVolume" = depth.simplicialVolume,
              "simplicial" = depth.simplicial,
              "Mahalanobis" = function(x, X) (.Mahalanobis_depth(x, colMeans(X), solve(cov(X)))),
              "projection" = depth.projection,
              "spatial" = depth.spatial
  )
  
  
  if (is.null(margins)){
    margins = c(min(data[,1]), max(data[,1]), min(data[,2]), max(data[,2]));  margins = margins + c(-0.1*(margins[2]-margins[1]), 0.1*(margins[2]-margins[1]), -0.1*(margins[4]-margins[3]), 0.1*(margins[4]-margins[3]))
  }
  
  gx <- seq(margins[1], margins[2], length=frequency)
  gy <- seq(margins[3], margins[4], length=frequency)
  y <- as.matrix(expand.grid(gx, gy))    
  
  depthcontours <- df(y,data,...)  
  
  contour(gx, gy, matrix(depthcontours, nrow=length(gx), ncol=length(gy)), add=TRUE, levels=depths*max(depthcontours), drawlabels=FALSE, col = col)
  
  return(depthcontours)
}


ellipsem <- function (mu, amat, c2, npoints = 100, showcentre = T, col, ...){
  if (all(dim(amat) == c(2, 2))) {
    eamat <- eigen(amat)
    hlen <- sqrt(c2/eamat$val)
    theta <- angle(eamat$vec[1, 1], eamat$vec[2, 1])
    ellipse(hlen[1], hlen[2], theta, mu[1], mu[2], npoints = npoints, col = col, ...)
    if (showcentre)
      points(mu[1], mu[2], pch = 3)
  }
  invisible()
}

ellipse <- function (hlaxa = 1, hlaxb = 1, theta = 0, xc = 0, yc = 0, newplot = F, npoints = 100, ...){
  a <- seq(0, 2 * pi, length = npoints + 1)
  x <- hlaxa * cos(a)
  y <- hlaxb * sin(a)
  alpha <- angle(x, y)
  rad <- sqrt(x^2 + y^2)
  xp <- rad * cos(alpha + theta) + as.double(xc)
  yp <- rad * sin(alpha + theta) + as.double(yc)
  if (newplot)
    plot(xp, yp, type = "l", ...)
  else lines(xp, yp, ...)
  invisible()
}

angle <- function (x, y)
{
  angle2 <- function(xy) {
    x <- xy[1]
    y <- xy[2]
    if (x > 0) {
      atan(y/x)
    }
    else {
      if (x < 0 & y != 0) {
        atan(y/x) + sign(y) * pi
      }
      else {
        if (x < 0 & y == 0) {
          pi
        }
        else {
          if (y != 0) {
            (sign(y) * pi)/2
          }
          else {
            NA
          }
        }
      }
    }
  }
  apply(cbind(x, y), 1, angle2)
}