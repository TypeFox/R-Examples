NULL
#
#' brickFromOutputSoil3DTensor
#' 
#' Extracts a brick or a raster layer from a output 3D Tensor or 2D map respectively 
#' 
#' @param x string. GEOtop keyword reletated to the 3D or 2D variable to be imported in R. 
#' @param when \code{\link{POSIXct-class}} for date and time on which the variable \code{x} is requested. 
#' @param layers number of soil layer or geotop keyword for soil layer (e.g. \code{SoilLayerThicknesses} or \code{SoilFile}). Default is  \code{SoilLayerThicknesses}. 
#' @param timestep time step expressed in seconds every which the raster file has been created. It can be a string corresponding to the geotop keyword in the inpts file. Default value is \code{"OutputSoilMaps"}.
#' @param suffix charachter string containing the decimal formatter used by GEOtop in the output file names. Default is "L\%04dN\%04.asc". A simple user is recommended not to modify the value of this argument and use the default value.
#' @param wpath,tz,use.read.raster.from.url see \code{\link{get.geotop.inpts.keyword.value}}
#' @param projfile  name of the \code{*.proj} file containing CRS information. See \code{\link{get.geotop.inpts.keyword.value}}. Default is \code{"geotop.proj"}. If is \code{NULL} or \code{NA} or this file does not exist, it is not searched and read.. In case \code{use.read.raster.from.url} is \code{TRUE} and no \code{NULL} or \code{NA} values are assinged, the \code{*.proj} file is searched. 
#' @param crs,start.from.zero see \code{\link{brick.decimal.formatter}}. If \code{crs} is not \code{NULL} (Default) , \code{projfile} is ignored.
#' @param one.layer logical value. If \code{TRUE} a \code{\link{RasterLayer-class}} object is imported, otherwise a \code{\link{RasterBrick-class}}object is returened. Default for \code{brickFromOutputSoil3DTensor} is \code{FALSE}
#' @param start_date_key,end_date_key initial and final detes and times of the GEOtop simulation or alternatively the respective keywords of \code{*.inpts} file (Default) 
#' @param secondary.suffix String secondary suffix which can be added at the end of the Map file name (optional). Default is \code{NULL} and no secondary suffix is added.  
#' @param ... additional arguments for \code{\link{get.geotop.inpts.keyword.value}} or \code{\link{brickFromOutputSoil3DTensor}}
#'  
#' @rdname brickFromOutputSoil3DTensor
#' @author Emanuele Cordano
#' 
#'  
#' 
#' @details These functions \code{brickFromOutputSoil3DTensor} and \code{rasterFromOutput2DMap} return 3D or 2D \code{\link{Raster-class}} objects respectively. \code{rasterFromOutput2DMap} is a wrapper function of \code{brickFromOutputSoil3DTensor} with the option \code{one.layer==TRUE}.
#' The functionswork with the following output keywords: 
#' 
#' \code{"SoilTempTensorFile"},
#' 
#' \code{"SoilAveragedTempTensorFile"},
#' 
#' \code{"SoilLiqContentTensorFile"},
#' 
#'\code{"SoilAveragedLiqContentTensorFile"},
#' 
#'\code{"SoilIceContentTensorFile"},
#' 
#' \code{"SoilAveragedIceContentTensorFile"},
#' 
#' \code{"SoilLiqWaterPressTensorFile"},
#' 
#'\code{"SoilTotWaterPressTensorFile"}     for \code{\link{brickFromOutputSoil3DTensor}};
#' 
#' 
#'\code{"FirstSoilLayerTempMapFile"}, 
#'   
#'\code{"FirstSoilLayerAveragedTempMapFile"},
#'  
#'\code{"FirstSoilLayerLiqContentMapFile"},
#' 
#'\code{"FirstSoilLayerIceContentMapFile"}, 
#'  
#'\code{"LandSurfaceWaterDepthMapFile"}, 
#'      
#'\code{"ChannelSurfaceWaterDepthMapFile"},
#'    
#'\code{"NetRadiationMapFile"}, 
#'               
#'\code{"InLongwaveRadiationMapFile"}, 
#'       
#'\code{"NetLongwaveRadiationMapFile"},
#'        
#'\code{"NetShortwaveRadiationMapFile"},
#'       
#'\code{"InShortwaveRadiationMapFile"}, 
#'       
#'\code{"DirectInShortwaveRadiationMapFile"},
#' 
#' 
#'\code{"ShadowFractionTimeMapFile"},
#' 
#'\code{"SurfaceHeatFluxMapFile"},  
#'           
#'\code{"SurfaceSensibleHeatFluxMapFile"},
#' 
#'\code{"SurfaceLatentHeatFluxMapFile"}, 
#'     
#'\code{"SurfaceTempMapFile"},    
#'      
#'\code{"PrecipitationMapFile"},  
#'             
#'\code{"CanopyInterceptedWaterMapFile"},
#'      
#'\code{"SnowDepthMapFile"},   
#'               
#'\code{"GlacierDepthMapFile"}, 
#'               
#'\code{"SnowMeltedMapFile"}, 
#'              
#' \code{"SnowSublMapFile"}, 
#'                   
#' \code{"GlacierMeltedMapFile"}, 
#'             
#' \code{"GlacierSublimatedMapFile"},  
#'         
#' \code{"AirTempMapFile"},        
#'             
#' \code{"WindSpeedMapFile"}, 
#' 
#' \code{"WindDirMapFile"}, 
#'                   
#' \code{"RelHumMapFile"}, 
#'    
#' \code{"SWEMapFile"},    
#'    
#' \code{"GlacierWaterEqMapFile"}  
#' 
#' \code{"SnowDurationMapFile"},
#'               
#' \code{"ThawedSoilDepthMapFile"},
#'             
#' \code{"ThawedSoilDepthFromAboveMapFile"}, 
#'   
#' \code{"WaterTableDepthMapFile"}, 
#' 
#' \code{"WaterTableDepthFromAboveMapFile"},  
#' 
#' \code{"NetPrecipitationMapFile"},
#' 
#' \code{"EvapotranspirationFromSoilMapFile"} for \code{\link{rasterFromOutput2DMap}}.
#' 

#' 
#' 
#' @seealso \code{\link{get.geotop.inpts.keyword.value}},\code{\link{brick.decimal.formatter}}
#'
#' @export
#' 
#' @examples 
#' 
#' library(geotopbricks)
#' # The data containing in the link are only for educational use
#' wpath <- "http://www.rendena100.eu/public/geotopbricks/simulations/idroclim_test1"
#' x <- "SoilLiqContentTensorFile"
#' tz <-  "Etc/GMT+1"
#' when <- as.POSIXct("2002-03-22",tz=tz)
#' 
#' # Not Run because it elapses too long time!!! 
#' # Please Uncomment the following lines to run by yourself!!!
#' # b <- brickFromOutputSoil3DTensor(x,when=when,wpath=wpath,tz=tz,use.read.raster.from.url=TRUE)
#' 
#' # a 2D map: 
#' x_e <- "SnowDepthMapFile"
#' # Not Run: uncomment the following line
#' # m <- rasterFromOutput2DMap(x_e,when=when,wpath=wpath,timestep="OutputSnowMaps",
#' #                            tz=tz,use.read.raster.from.url=TRUE)
#' ## NOTE: set use.read.raster.from.url=FALSE (default) 
#' # if the "wpath" directorty is in the local file system.
#' # Not Run: uncomment the following line
#' # plot(m)
#'

#
#SoilLiqContentTensorFile
#
#
#SoilLiqContentTensorFile 
#
#Name of the ensamble of raster maps corresponding to the liquid water
#content of each layer (if PlotSoilDepth6=0 it writes the value at the
#corresponding depths)
#
#SoilLiqWaterPressTensorFile 
#
#
#
#SoilTempTensorFile Name of the ensamble of raster maps corresponding to the temperature
#of each layer (if PlotSoilDepth6=0 it writes the value at the corresponding
#depths)
#SoilAveragedTempTensorFile Name of the ensamble of raster maps corresponding to the average temperature
#of each layer (if PlotSoilDepth6=0 it writes the value at the
#corresponding depths)
#IceLiqContentTensorFile Name of the ensamble of raster maps corresponding to the average ice
#content of each layer (if PlotSoilDepth6=0 it writes the value at the corresponding
#depths)
#Table 15.10: Keywords of print output tensor maps for soil and
#
#
#
#
##




brickFromOutputSoil3DTensor <- function(x,when,layers="SoilLayerThicknesses",one.layer=FALSE,suffix="L%04dN%04d.asc",wpath=NULL,tz="A",start_date_key="InitDateDDMMYYYYhhmm",end_date_key="EndDateDDMMYYYYhhmm",timestep="OutputSoilMaps",use.read.raster.from.url=FALSE,crs=NULL,projfile="geotop.proj",start.from.zero=FALSE,secondary.suffix=NULL,...) {
	
	out <- NULL

	
	if (is.null(crs)) { 
		
		cond <- (!is.null(wpath)) | is.null(projfile) | is.na(projfile)
		if (cond) { 
			projfile <- paste(wpath,projfile,sep="/")
			
			
			cond <- file.exists(projfile) | use.read.raster.from.url
			crs <- getProjection(projfile,cond=cond)
		#	if (cond) {
		#	#	crs <- readLines(projfile,warn=FALSE)
		#		crs <- scan(projfile,what="list",sep="\n",n=1)
		#		
		#	}
		}
	}
	
	
	
	
	start_s <- get.geotop.inpts.keyword.value(start_date_key,date=TRUE,wpath=wpath,tz=tz,...) ###wpath=wpath,tz="A")
	end_s <- get.geotop.inpts.keyword.value(end_date_key,date=TRUE,wpath=wpath,tz=tz,...) ###wpath=wpath,tz="A")
	
	
	if (!is.numeric(timestep)) timestep <- get.geotop.inpts.keyword.value(timestep,wpath=wpath,numeric=TRUE,...)*3600 
	
	time <- seq(from=start_s,to=end_s,by=timestep)
	
	
	
	if (is.numeric(layers)) {
	
		if ((length(layers)==1) & (layers[1]>1)) layers <- 1:layers 
		if ((length(layers)==1) & (layers[1]=1)) one.layer <- TRUE
	
	} else {
		
		if (layers=="SoilLayerThicknesses") {
		
			###get.geotop.inpts.keyword.value("SoilLayerThicknesses",numeric=TRUE,wpath=wpath,...)
			layers <- get.geotop.inpts.keyword.value("SoilLayerThicknesses",numeric=TRUE,wpath=wpath,...) ####wpath=wpath)
			###print(layers)
			if (is.null(layers)) {
				layers <-  "SoilParFile"
			} else {
				
				#layers <- 1:length(layers)
			}	
				
		} 
		
		if (layers[1]=="SoilParFile") {
			
			layers <- get.geotop.inpts.keyword.value("SoilParFile",wpath=wpath,add_wpath=TRUE,data.frame=TRUE,level=1,date_field=NULL,...)$Dz
			
		
			

		} else if (!is.numeric(layers)) {
			
			layers <- 1:2
			
		}
		
		
	}
	
	map.prefix <- get.geotop.inpts.keyword.value(x,numeric=FALSE,date=FALSE,wpath=wpath,add_wpath=TRUE,...)
	
	
	if (!is.null(secondary.suffix)) {
			
		map.prefix <- paste(map.prefix,secondary.suffix,sep="")
		
	}
	
	
	
	
	
	print(paste("Maps to import:",length(when),"from",as.character(when[1]),"to",as.character(when[length(when)]),sep=" "))
	
	out <- lapply(X=when,FUN=function(whenx,map.prefix,suffix,crs,layers,start.from.zero,one.layer,time,timestep) {
		
		print(paste("Importing",as.character(whenx),sep=" "))
		t_index <- abs(whenx-time)<timestep
		index <- 1
		n <- which(t_index)[index]
				
		if (length(n)>0) {
					
				n <- n[1]
					
		} 		
				
				
		if (!is.na(n)) {
		
			suffix <- str_split(suffix,"N")[[1]]
			suffix[2] <- sprintf(suffix[2],n)
		
			if (one.layer) {
		
				suffix <- paste("N",suffix[2],sep="")
			
			
			} else {	
				suffix <- paste(suffix,collapse="N")
			}	
		
			map.filename <- paste(map.prefix,suffix,sep="")
		
			if (one.layer) {
			
				if (use.read.raster.from.url) {				
					out <- read.raster.from.url(x=map.filename)
				} else {
					out <- raster(map.filename)
				}
			
				if (!is.null(crs)) {
					projection(out) <- crs
				}
			
			
			} else {
			
				out <- brick.decimal.formatter(map.filename,nlayers=length(layers),use.read.raster.from.url=use.read.raster.from.url,start.from.zero=start.from.zero,crs=crs)
		
				if (start.from.zero) {
			
					names(out) <- paste("L",0:length(layers),sep="")
				} else {
					names(out) <- paste("L",1:length(layers),sep="")
			
				}
		
				names(out)[1:length(layers)] <- paste(names(out)[1:length(layers)],layers,sep="_")
			} 

		
		} else {
		
		out <- NULL
		}
	
		return(out)
	},map.prefix=map.prefix,suffix=suffix,crs=crs,layers=layers,start.from.zero=start.from.zero,one.layer=one.layer,time=time,timestep=timestep)
	
	
	
	names(out) <- strftime(when,tz=tz,format="DATE-TIME %Y-%m-%d %H:%M")
	
	
	if (length(out)==1) out <- out[[1]]
	
	return(out)
	
	
}

NULL 
#'
#' 
#' @rdname brickFromOutputSoil3DTensor
#' @export 

rasterFromOutput2DMap <- function(x,when,...) {
	
	
	out <- brickFromOutputSoil3DTensor(x=x,when=when,layers=1,one.layer=TRUE,...)
	return(out)
	
}


