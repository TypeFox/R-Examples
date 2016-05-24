#'LiDAR-derived Canopy Height Model (CHM) smoothing
#'
#'@description LiDAR-derived Canopy Height Model (CHM) smoothing is used to eliminate spurious local maxima caused by tree branches.
#'
#'@usage CHMsmoothing(chm, filter, ws, sigma)
#'
#'@param chm A LiDAR-derived Canopy Height Model (CHM) RasterLayer or SpatialGridDataFrame file.
#'@param ws The dimension of a window size, e.g. 3,5, 7 and so on. Default is 5.
#'@param filter Filter type: mean, median, maximum or Gaussian. Default is mean.
#'@param sigma Used only when filter parameter is equal to Gaussian, e.g. 0.5, 1.0, 1.5 and so on. Default is 0.67. 
#'@return Returns a CHM-smoothed raster.
#'@author Carlos Alberto Silva. 
#'@seealso \code{\link[raster]{focal}} in the \emph{raster} package.
#'@examples
#'
#'#=======================================================================#
#'# Importing the LiDAR-derived CHM file
#'data(chm) # or set a CHM. e.g. chm<-raster("CHM_stand.asc") 
#'
#'#=======================================================================#
#'# Example 01: Smoothing the CHM using a Gaussian filter
#'#=======================================================================#
#'# Set the ws:
#'ws<-3 # dimension 3x3
#'
#'# Set the filter type
#'filter<-"Gaussian"
#'
#'# Set the sigma value
#'sigma<-0.6
#'
#'# Smoothing CHM
#'sCHM<-CHMsmoothing(chm, filter, ws, sigma)
#'
#'#=======================================================================# 
#'# Example 02: Smoothing the CHM using a mean filter
#'#=======================================================================#
#'# Set the ws:
#'ws<-5 # dimension 5x5
#'
#'# Set the filter type
#'filter<-"mean"
#'
#'# Smoothing and plotting LiDAR-derived CHM 
#'sCHM<-CHMsmoothing(chm, filter, ws)
#'
#'@importFrom raster raster focal
#'@export CHMsmoothing
CHMsmoothing<-function(chm, filter="mean", ws=5, sigma=0.6) {

  if (class(chm)[1]!='RasterLayer') {
      chmInput<-as(chm,'RasterLayer')
      } else {chmInput<-chm
  }
  
  if (filter!="mean" & filter!="median" & filter!="maximum" & filter!="Gaussian") {stop("The filter parameter is invalid. Please, use one of this filter types: 'mean','median','maximum','minimum',gaussian'")}
  if (class(ws)!="numeric") {stop("The ws parameter is invalid. It is not a numeric input")}
  if (class(sigma)!="numeric") {stop("The sigma parameter is invalid. It is not a numeric input")}
  
  if (filter == "mean") {
    wf<-matrix(c(rep(1,ws*ws)),nrow=ws,ncol=ws)
    chmR <- focal(chmInput, w=wf, fun=mean)
  }
  if (filter == "median") {
    wf<-matrix(c(rep(1,ws*ws)),nrow=ws,ncol=ws)
    chmR <- focal(chmInput, w=wf, fun=median)
  }
  if (filter == "maximum") {
    wf<-matrix(c(rep(1,ws*ws)),nrow=ws,ncol=ws)
    chmR <- focal(chmInput, w=wf, fun=max)
  }
  
  if (filter =="Gaussian") {
    
    fgauss <- function(sigma, n=ws) {
      m <- matrix(ncol=n, nrow=n)
      col <- rep(1:n, n)
      row <- rep(1:n, each=n)
      x <- col - ceiling(n/2)
      y <- row - ceiling(n/2)
      m[cbind(row, col)] <- 1/(2*pi*sigma^2) * exp(-(x^2+y^2)/(2*sigma^2))
      m / sum(m)
    }
    gf=fgauss(sigma)
   
    chmR <- focal(chmInput, w=gf)
  }
  
  return(chmR)
  
}
