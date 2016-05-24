###' lgcp
###' 
###' An R package for spatiotemporal prediction and forecasting for log-Gaussian Cox processes.
###' 
###' \tabular{ll}{
###' Package: \tab lgcp\cr
###' Version: \tab 0.9-4-1\cr
###' Date: \tab 2012-17-04\cr
###' License: \tab GPL (>= 2) \cr
###' }
###' 
###'  
###' For examples and further details of the package, type vignette("lgcp"), or refer to the paper associated with this package.
###' 
###' The content of \code{lgcp} can be broken up as follows:\cr
###' 
###' \emph{Datasets} wpopdata.rda, wtowncoords.rda, wtowns.rda. Giving regional and town poopulations as well as town coordinates,are provided by Wikipedia 
###' and The Office for National Statistics under respectively
###' the Creative Commons Attribution-ShareAlike 3.0 Unported License and the Open Government Licence.\cr
###' 
###' \emph{Data manipulation}\cr
###' 
###' \emph{Model fitting and parameter estimation}\cr
###' 
###' \emph{Unconditional and conditional simulation}\cr
###' 
###' \emph{Summary statistics, diagnostics and visualisation}\cr
###' 
###' }
###' 
###' \section{Dependencies}{
###' The \code{lgcp} package depends upon some other important contributions to CRAN in order to operate; their uses here are indicated:\cr\cr
###'     spatstat, sp, RandomFields, iterators, ncdf, methods, tcltk, rgl, rpanel, fields, rgdal, maptools, rgeos, raster
###' }
###' 
###' \section{Citation}{
###' To see how to cite \code{lgcp}, type \code{citation("lgcp")} at the console.
###' }
###' 
###' \references{
###' \enumerate{
###'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
###'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
###'     \item Wood ATA, Chan G (1994). Simulation of Stationary Gaussian Processes in [0,1]d. Journal of Computational and Graphical Statistics, 3(4), 409-432.
###'     \item Moller J, Syversveen AR, Waagepetersen RP (1998). Log Gaussian Cox Processes. Scandinavian Journal of Statistics, 25(3), 451-482.
###' }
###' 
###' @docType package
###' @name lgcp-package
###' @author Benjamin Taylor, Health and Medicine, Lancaster University,
###'  Tilman Davies, Institute of Fundamental Sciences - Statistics, Massey University, New Zealand., 
###'  Barry Rowlingson, Health and Medicine, Lancaster University
###'  Peter Diggle, Health and Medicine, Lancaster University
###' @keywords package


###' @import methods
###' @importFrom RandomFields CovarianceFct
###' @importFrom Matrix sparseMatrix
###' @importFrom rpanel rp.block rp.button rp.control rp.radiogroup rp.slider rp.textentry rp.checkbox
###' @importFrom tcltk setTkProgressBar tclvalue tkProgressBar tkwinfo tkwm.geometry 
###' @importFrom rgeos gDisjoint gIntersection gIntersects gTouches gUnaryUnion 
###' @importFrom ncdf4 nc_open nc_close nc_sync ncvar_get ncdim_def ncvar_def nc_create ncvar_put
###' @importFrom maptools label
###' @importFrom iterators icount iter nextElem 
###' @importFrom fields image.plot  
###' @importFrom tcltk setTkProgressBar tclvalue tkProgressBar tkwinfo tkwm.geometry 
###' @importFrom raster raster disaggregate aggregate resample brick as.data.frame
###' @importFrom sp proj4string<- proj4string SpatialPixelsDataFrame SpatialGridDataFrame Polygon Polygons SpatialPolygons coordinates CRS geometry GridTopology over proj4string SpatialGrid SpatialPixels SpatialPoints SpatialPolygonsDataFrame split spTransform 
###' @importFrom spatstat rescale.ppp resolve.2D.kernel convexhull density.ppp setminus.owin simplify.owin interp.im affine.ppp rescale as.ppp as.owin affine area area.owin as.im as.mask as.polygonal as.rectangle im inside.owin is.im is.multitype is.polygonal is.ppp is.rectangle Kest Kinhom marks nearest.raster.point owin pairdist pcf pcfinhom ppp rmax.rule spatstat.options termsinformula variablesinformula verifyclass 
###' @importFrom stats acf approx as.formula C coefficients density deriv df dlnorm dnorm end fft formula Gamma glm lm lowess median model.matrix optim optimise quantile rnorm rpois runif sd start terms vcov window quasipoisson
###' @importFrom graphics abline grid hist identify image legend lines locator par plot points title
###' @importFrom grDevices dev.off devAskNewPage heat.colors rainbow rgb
###' @importFrom utils assignInNamespace browseURL flush.console head menu object.size setTxtProgressBar tail txtProgressBar

`lgcp` = NA

# package ncdf obsolete
### @importFrom ncdf close.ncdf get.var.ncdf open.ncdf dim.def.ncdf var.def.ncdf create.ncdf put.var.ncdf sync.ncdf