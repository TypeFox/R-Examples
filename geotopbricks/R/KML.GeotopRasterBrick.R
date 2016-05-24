NULL
#' \code{KML} method for a \code{GeotopRasterBrick} object
#' 
#' @param x the \code{\link{GeotopRasterBrick}} object
#' @param filename mane of the KML file to produce
#' @param crs character string containg the LatLon reference system. Default is \code{"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"} (see \url{http://spatialreference.org/ref/epsg/4326/}). 
#' @param ... further argument for S4 method \code{KLM} for Raster object.
#' 
#' @note A coordinate transformation is made with \code{\link{projectRaster}}. 
#' 
#' @title KML
#' @name KML
#' 
#' @export
#' @rdname KML-methods
#' @keywords methods
#' @docType methods
#' @method KML GeotopRasterBrick 
#' @aliases KML,GeotopRasterBrick-method

#' @examples
#' 
#' 
#' library(geotopbricks)
#' # The examples is the following R script conteined in a 'inst' directory of the package source
#' f <- system.file("doc/examples/example.KML.GeotopRasterBrick.R",package="geotopbricks")
#' #  source(f) # Uncomment this line to run the example. 
#' # You can copy the example file using file.copy(from=f,to=....,...) See file.copy documentation


### NO \code{as.charachter("+init=epsg:4326")}(see URL)
# rm(list=ls())
# library(rgdal)
# library(raster)
# library(zoo)
# library(geotopbricks)
#
# ## working path for 3D distributed raster maps 
# ## The study case of Ton-Toss (Val di Non, Trentino, Italy)
#
# wpath <- 'http://meteogis.fmach.it/idroclima/ton-toss' 
# # WARNING: In order to save disk space, some files of this simulation (unusuful for the example) were removed !!!!
# # keyword for water content 3D+time output raster maps
# watercontent_prefix <- get.geotop.inpts.keyword.value("SoilLiqContentTensorFile",wpath=wpath) #"thetaliq"
# #  crs projection
# crs  <-"+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "
#
# #  vector with vertical layer thckness  
#
# layers <- get.geotop.inpts.keyword.value("SoilLayerThicknesses",numeric=TRUE,wpath=wpath) 
# names(layers) <- paste("L",1:length(layers))
#
# # set time during wich GEEOtop simulation provided maps (the script is written for daily frequency")
#
# start <-  get.geotop.inpts.keyword.value("InitDateDDMMYYYYhhmm",date=TRUE,wpath=wpath,tz="A")
# end <- get.geotop.inpts.keyword.value("EndDateDDMMYYYYhhmm",date=TRUE,wpath=wpath,tz="A")
#
# # set time during wich GEEOtop simulation provided maps (the script is written for daily frequency")
#
#
# # In this examples maps are provided with daily frequency!!
# time <- seq(from=start,to=end,by="days") 
#
# #
# # files extracts the filename of the maps of the first 4 layers!! 
# # 
# files <- pointer.to.maps.xyz.time(wpath=wpath,map.prefix=watercontent_prefix,zoo.index=time,nlayers=4)
#
# # import maps 
#
# start_m <- as.POSIXlt("2012-04-15 00:00",tz="A")
# end_m <- as.POSIXlt("2012-04-20 00:00",tz="A")
#
# gt_wtc1 <- geotopbrick(x=files,layer=1,timerange=c(start_m,end_m),crs=crs)
# gt_wtc2 <- geotopbrick(x=files,layer=2,timerange=c(start_m,end_m),crs=crs)
# gt_wtc3 <- geotopbrick(x=files,layer=3,timerange=c(start_m,end_m),crs=crs)
#
# # Averaged Soil Water Content in the first 3 soil layers (about 33 centimetes)
# # 'RasterBrick objects'
# wtc <- (brick(gt_wtc1)*layers[1]+brick(gt_wtc2)*layers[2]+brick(gt_wtc3)*layers[3])/sum(layers[1:3])
# # 'GeotopRasterBrick objects
# gt_wtc <- ((gt_wtc1)*layers[1]+(gt_wtc2)*layers[2]+(gt_wtc3)*layers[3])/sum(layers[1:3])
#
# # set colors for 'KML' and 'plot' 
#
# N <- 10000
#
# start_watercontent <- 4/12 
# end_watercontent <- 9/12
# col_watercontent <- rainbow(start=start_watercontent,end=end_watercontent,n=N,alpha=0.6)
# # Creates the KML 
# KML(gt_wtc,filename="zz_wtc_33cm_ton_toss.kml",overwrite=TRUE,col=col_watercontent)	
#
# # raster legend (in 
# pdf <- "zz_wtc_33cm_ton_toss_legend.pdf"
# color.bar.raster(x=gt_wtc,col_watercontent,digits=2,pdf=pdf)
#
#






# wc <- projectRaster(watercontent_a,crs="+init=epsg:4326")
#> raster2 <- projectRaster(raster1,CRS("+init=epsg:4326"))
#plot
#Error in function (classes, fdef, mtable)  : 
#			unable to find an inherited method for function "res", for signature "CRS"
#In addition: Warning message:
#		In min(dim(to)[1:2]) : no non-missing arguments to min; returning Inf
#> raster2 <- projectRaster(raster1,crs=CRS("+init=epsg:4326"))
#Error in .computeRes(from, projto) : 
#		STRING_ELT() can only be applied to a 'character vector', not a 'S4'
#> raster2 <- projectRaster(raster1,crs="+init=epsg:4326"))
#Error: unexpected ')' in "raster2 <- projectRaster(raster1,crs="+init=epsg:4326"))"
#> raster2 <- projectRaster(raster1,crs="+init=epsg:4326")
#> KML(raster2,filename="prova.kml")
#> wc <- projectRaster(watercontent_a,crs="+init=epsg:4326")
#> KML(wc,filename="provawc.kml")
#> wc <- projectRaster(watercontent_a,crs="+init=epsg:4326")

setMethod('KML', signature(x='GeotopRasterBrick'), 
		function (x, filename, crs=as.character("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"),...) {  # "+init=epsg:4326")#} ,zip='', overwrite=FALSE, ...) {
		
			y <- projectRaster(brick(x),crs=crs)			
			out <- KML(x=y,filename=filename,...) 
			return(out)
		}

)
