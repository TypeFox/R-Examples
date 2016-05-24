NULL
#'
#' \code{plot} method for a \code{GeotopRasterBrick} object
#' 
#' @param x the \code{\link{GeotopRasterBrick}} object
#' @param y further argument
#' @param ... further argument for S4 method \code{plot} for Raster object.
#' 
#' 
#' @title plot
#' @name plot
#' @export
#' @rdname plot-methods
#' @keywords methods
#' @docType methods
#' @method plot GeotopRasterBrick 
#' @aliases plot,GeotopRasterBrick,ANY-method
#' 

#' 
#' @seealso \code{\link{KML}}
#' 
#' @examples 
#' 
#' 
#' library(geotopbricks)
#' # The examples is the following R script conteined in a 'inst' directory of the package source
#' f <- system.file("doc/examples/example.plot.GeotopRasterBrick.R",package="geotopbricks")
#' #  source(f) # Uncomment this line to run the example. 
#' # You can copy the example file using file.copy(from=f,to=....,...) See file.copy documentation

#rm(list=ls())
#library(rgdal)
#library(raster)
#library(zoo)
#library(geotopbricks)
#'
### working path for 3D distributed raster maps 
### The study case of Ton-Toss (Val di Non, Trentino, Italy)
#
#wpath <- 'http://meteogis.fmach.it/idroclima/ton-toss' 
## WARNING: In order to save disk space, some files of this simulation (unuseful for the example) were removed !!!!
## keyword for water content 3D+time output raster maps
#watercontent_prefix <- get.geotop.inpts.keyword.value("SoilLiqContentTensorFile",wpath=wpath) #"thetaliq"
##  crs projection
#crs  <-"+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "
#
##  vector with vertical layer thckness  
#
#layers <- get.geotop.inpts.keyword.value("SoilLayerThicknesses",numeric=TRUE,wpath=wpath) 
#names(layers) <- paste("L",1:length(layers))
#'
## set time during wich GEEOtop simulation provided maps (the script is written for daily frequency")
#'
#start <-  get.geotop.inpts.keyword.value("InitDateDDMMYYYYhhmm",date=TRUE,wpath=wpath,tz="A")
#end <- get.geotop.inpts.keyword.value("EndDateDDMMYYYYhhmm",date=TRUE,wpath=wpath,tz="A")
#'
## set time during wich GEEOtop simulation provided maps (the script is written for daily frequency")
#
#
## In this examples maps are provided with daily frequency!!
#time <- seq(from=start,to=end,by="days") 
#
#
## files extracts the filename of the maps of the first 4 layers!! 
## 
#files <- pointer.to.maps.xyz.time(wpath=wpath,map.prefix=watercontent_prefix,zoo.index=time,nlayers=4)
#
## import maps 
#
#start_m <- as.POSIXlt("2012-04-15 00:00",tz="A")
#end_m <- as.POSIXlt("2012-04-20 00:00",tz="A")
#
#gt_wtc1 <- geotopbrick(x=files,layer=1,timerange=c(start_m,end_m),crs=crs)
#gt_wtc2 <- geotopbrick(x=files,layer=2,timerange=c(start_m,end_m),crs=crs)
#gt_wtc3 <- geotopbrick(x=files,layer=3,timerange=c(start_m,end_m),crs=crs)
#'
## Averaged Soil Water Content in the first 3 soil layers (about 33 centimetes)
## 'RasterBrick objects'
#wtc <- (brick(gt_wtc1)*layers[1]+brick(gt_wtc2)*layers[2]+brick(gt_wtc3)*layers[3])/sum(layers[1:3])
## 'GeotopRasterBrick objects
#gt_wtc <- ((gt_wtc1)*layers[1]+(gt_wtc2)*layers[2]+(gt_wtc3)*layers[3])/sum(layers[1:3])
#'
## set colors for 'KML' and 'plot' 
#'
#N <- 10000
#
#start_watercontent <- 4/12 
#end_watercontent <- 9/12
#col_watercontent <- rainbow(start=start_watercontent,end=end_watercontent,n=N,alpha=1.0)
#
#  
### plot(gt_wtc,col=col_watercontent) - Do Not Run - Uncomment to run the lines!
#

#

setMethod('plot', signature(x='GeotopRasterBrick',y='ANY'), 
		function (x,y=NULL,...) { #} ,zip='', overwrite=FALSE, ...) {
		
			
			out <- plot(brick(x),...) 
			return(out)
		}

)
