##' spatsurv
##' 
##' An R package for spatially correlated parametric proportional hazards survial analysis. 
##' 
##' \tabular{ll}{
##' Package: \tab spatsurv\cr
##' Version: \tab 0.9-11\cr
##' Date: \tab 2015-06-24\cr
##' License: \tab GPL-3 \cr
##' }
##' 
##'  
##' 
##' section{Dependencies}{
##' The package \code{spatsurv} depends upon some other important contributions to CRAN in order to operate; their uses here are indicated:\cr\cr
##'     survival, sp, spatstat, raster, iterators, RandomFields, fields, rgl, Matrix, stringr, RColorBrewer, geostatsp, rgeos.
##' }
##' 
##' section{Citation}{
##' To cite use of \code{spatsurv}, the user may refer to the following work:\cr\cr
##' spatsurv: an \code{R} Package for Bayesian Inference with Spatial Survival Models.\cr 
##' Benjamin M. Taylor and Barry S. Rowlingson. \cr
##' Submitted to The Journal Of Statistical Software.
##' }
##' 
##' references{
##' X
##' }
##' 
##' @docType package
##' @name spatsurv-package
##' @author Benjamin Taylor, Health and Medicine, Lancaster University,
##'  Barry Rowlingson, Health and Medicine, Lancaster University
##' @keywords package
##'
##'



##' @importFrom OpenStreetMap openmap
##' @importFrom RColorBrewer brewer.pal 
##' @importFrom stringr str_count str_detect
##' @importFrom Matrix Matrix sparseMatrix
##' @importFrom rgl abclines3d aspect3d axes3d planes3d points3d segments3d text3d title3d 
##' @importFrom fields image.plot  
##' @importFrom RandomFields CovarianceFct
##' @importFrom rgeos gBuffer
##' @importFrom iterators icount iter nextElem
##' @importFrom sp bbox proj4string<- proj4string SpatialPixelsDataFrame SpatialGridDataFrame Polygon Polygons SpatialPolygons coordinates CRS geometry GridTopology over proj4string SpatialGrid SpatialPixels SpatialPoints SpatialPolygonsDataFrame split spTransform 
##' @importFrom spatstat rpoint progressreport
##' @importFrom survival Surv survfit
##' @importFrom geostatsp asImRaster
##' @importFrom raster crop brick raster


## @import stats
## @import graphics 
## @import methods
## @import utils
## @import grDevices

##' @importFrom stats as.formula acf coefficients deriv dexp dist dnorm end fft fitted formula Gamma integrate knots lm model.matrix optim optimise poly quantile rbinom rexp rnorm runif sd start update var residuals cov
##' @importFrom graphics hist legend lines matplot par plot points title abline
##' @importFrom methods as
##' @importFrom utils txtProgressBar setTxtProgressBar browseURL flush.console
##' @importFrom grDevices adjustcolor 





`spatsurv` = NA


