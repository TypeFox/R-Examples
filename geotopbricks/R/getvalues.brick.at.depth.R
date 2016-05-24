# TODO: Add comment
# 
# Author: ecor
###############################################################################
NULL
#'
#' Interpolates the values of a 'brick' at a certain depth  and returns the map of brick values at the "depth" level
#' 
#' @param x a 'RasterBrick' or a three-dimensional array
#' @param depth depth map, generally a 'RasterLayer' object
#' @param layers vector of layer thickness
#' @param i0 a 'Raster' containing the number of soil laver just over the bedrock. Default is \code{NULL} and is then calculated.
#' @param verify logical. Default is \code{FALSE}. If it is \code{TRUE}, it verifies that function is working correctly. 
#' @param ... further argument 
#' 
#' @note \code{x} and \code{depth} or \code{i0} must cover the same spatial region.
#' 
#' @return a list of 'Raster' maps: 
#' 
#' \code{i0}   a 'Raster' containing the number of soil laver just over the bedrock
#' 
#' \code{val_z0}   a 'Raster' containing the values of \code{x} at the \code{i0}-th layer
#'
#' \code{val_z1}   a 'Raster' containing the values of \code{x} at the (\code{i0}+1)-th layer 
#' 
#' \code{z0}       a 'Raster' containing the depth of the center  of  the \code{i0}-th layer
#' 
#' \code{z1}       a 'Raster' containing the depth of the center  of  the (\code{i0}+1)-th layer
#' 
#' 
#' @seealso  code{\link{vertical.aggregate.brick.within.depth}}
#' 
#' @export
#' 
#' @examples 
#' 
#' library(geotopbricks)
#' # The examples is the following R script conteined in a 'inst' directory of the package source
#' f <- system.file("doc/examples/example.getvalues.brick.at.depth.R",package="geotopbricks")
#' #  source(f) # Uncomment this line to run the example. 
#' # You can copy the example file using file.copy(from=f,to=....,...) See file.copy documentation
# library(rgdal)
# library(raster)
# library(zoo)
# library(methods)
# library(geotopbricks)
# wpath <- 'http://meteogis.fmach.it/idroclima/panola13_run2xC_test3'  
# 
# prefix <- get.geotop.inpts.keyword.value("SoilLiqWaterPressTensorFile",wpath=wpath)
#
# 
# slope <- get.geotop.inpts.keyword.value("SlopeMapFile",raster=TRUE,wpath=wpath) 
# bedrock_depth <- get.geotop.inpts.keyword.value("BedrockDepthMapFile",raster=TRUE,wpath=wpath) 
# 
# 
# 
# layers <- get.geotop.inpts.keyword.value("SoilLayerThicknesses",numeric=TRUE,wpath=wpath)
# 
# names(layers) <- paste("L",1:length(layers),sep="")
# 
# # set van genuchten parameters to estimate water volume 
# theta_sat <- get.geotop.inpts.keyword.value("ThetaSat",numeric=TRUE,wpath=wpath)
# theta_res <- get.geotop.inpts.keyword.value("ThetaRes",numeric=TRUE,wpath=wpath)
# alphaVG <-  get.geotop.inpts.keyword.value("AlphaVanGenuchten",numeric=TRUE,wpath=wpath) # expressed in mm^-1
# nVG <-  get.geotop.inpts.keyword.value("NVanGenuchten",numeric=TRUE,wpath=wpath) 
# 
# 
# # end set van genuchten parameters to estimate water volume
#
# 
# # set time during wich GEEOtop simulation provided maps (the script is written for daily frequency")
# 
# start <-  get.geotop.inpts.keyword.value("InitDateDDMMYYYYhhmm",date=TRUE,wpath=wpath,tz="A") 
# end <- get.geotop.inpts.keyword.value("EndDateDDMMYYYYhhmm",date=TRUE,wpath=wpath,tz="A") 
# 
# # end set time during wich GEEOtop simulation provided maps (the script is written for daily frequency")
# 
# ## The maps are obtanied with daily frequancy
#
# time <- seq(from=start,to=end,by="days") 
#
# ## Pressure head map filename
#
# psiliq_files <- pointer.to.maps.xyz.time(wpath=wpath,map.prefix=prefix,zoo.index=time,nlayers=length(layers))
#
# ## Note that in this similation 'psi' maps are returned with daily frequency!!!!
# days <- c(1,2,5,10,15,20) # integer!!! 
#
# i <- 5 # ten days since the beginning of the simulation period 
# 
# psi <- brick(psiliq_files,layer=1:length(layers),rows=days[i])
# psi_bedrock <- getvalues.brick.at.depth(psi,depth=bedrock_depth,layers=layers,verify=FALSE)
# # Plot the map of the presuure head over the bedrock:
# # plot(psi_bedrock$val_z0) # expressed in millimeters !!!! - Do not Run
# 
# ## Calculation of water-table thickness according to Cordano and Rigon's Equation 10 (http://www.agu.org/pubs/crossref/2008/2006WR005740.shtml) 
# 
# slope_cos <- cos(slope*pi/180)
# 
# ## geotop_watertable_thickness: i.e. the slope-normal thickness of water-table film over the bedrock!!! 
# 
# geotop_watertable_thickness <- (psi_bedrock$val_z0/slope_cos+bedrock_depth-psi_bedrock$z0)
# geotop_watertable_thickness[geotop_watertable_thickness<0] <- 0
# 
# # plot(geotop_watertable_thickness) # expressed in millimeters !!!! - Do not Run
# 
#







getvalues.brick.at.depth<- function (x,depth,layers,i0=NULL,verify=FALSE,...) {
	

	if (nlayers(x)!=length(layers)) {
		
		print("Error in getvalues.brick.at.depth: dimensions mismatching between x and layers ")
		return(-1)
	}
	
	if (nrow(x)!=nrow(depth)) {
		
		print("Error in getvalues.brick.at.depth: dimensions (rows) mismatching between x and depth ")
		return(-1)
	}
	
	if (ncol(x)!=ncol(depth)) {
		
		print("Error in getvalues.brick.at.depth: dimensions (columns) mismatching between x and depth ")
		return(-1)
	}
	

	
	z <- layers/2
	
	
	
	
	for (i in 2:length(z)) {
		
		z[i] <- z[i-1]+(layers[i-1]+layers[i])/2
	}
	
	
	if (is.null(i0)) {
		
		i0 <- depth*0+0

	
	
		for (i in 1:length(z)) {
		#	print(i)
			i0[depth>z[i]] <- i0[depth>z[i]]+1
		}

	
	#	i0[i0<=1] <- 1 
	#	i0[i0>length(z)-1] <- length(z)-1


	} else if (nrow(x)!=nrow(i0)) {
		
		print("Error in getvalues.brick.at.depth: dimensions (rows) mismatching between x and i0 ")
		return(-1)
	} else if (ncol(x)!=ncol(i0)) {

		print("Error in getvalues.brick.at.depth: dimensions (columns) mismatching between x and i0 ")
		return(-1)
	} 
	
		
		
	i0[i0<=1] <- 1 
	i0[i0>length(z)-1] <- length(z)-1	
	
	out <- list()
	
	out$i0 <- i0
	out$val_z0 <- depth*0.0
	out$val_z1<- depth*0.0
	out$z0 <- depth*0.0 
	out$z1 <- depth*0.0
	
	for (i in 1:(length(z)-1)) {
	
		temp0 <- subset(x,i)
		temp1 <- subset(x,i+1)
		out$z0[out$i0==i] <- z[i]
		out$z1[out$i0==i] <- z[i+1]
		out$val_z0[out$i0==i] <- temp0[out$i0==i]
		out$val_z1[out$i0==i] <- temp1[out$i0==i]
		
	}	
	
	
	
#	
#	for (r in 1:nrow(i0)) {
#		for (c in 1:ncol(i0)) {
#	
#			kindex <- i0[r,c]
#			if (!is.na(kindex)) {
#			
#			print
#				out$val_z0[r,c] <- x[kindex,r,c]
#				out$val_z1[r,c] <- x[kindex+1,r,c]
#				out$z0[r,c] <- z[kindex] 
#				out$z1[r,c]  <- z[kindex+1]
#			}
#		
#		}
#		
#	}

	if (verify) {
		
		zmin <- min(z,na.rm=TRUE)
		zmax <- max(z,na.rm=TRUE)
		imax <- length(z)-2
		imin <- 1 

		verified <- as.matrix((depth>=out$z0) & (depth<=out$z1))
		verified <- as.matrix((depth>=zmax) & (out$i0==imax)) | verified
		verified <- as.matrix((depth<=zmin) & (out$i0==imin)) | verified
	
		xverified <- min(verified,na.rm=TRUE)
	    if (xverified!=1)	{
			
			 print("Error in getvalues.brick.at.depth: depth bedrock is not included between z0 and z1, function does not work correctly!")
			 return(-1)
			 
		 }
	}
	return(out)
}


