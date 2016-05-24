# TODO: Add comment
# 
# Author: ecor
###############################################################################
NULL
#'
#' Aggregates with a mean or an addition on the vertical profile the values of a 'brick' within a certain depth  and returns the vertical aggregated map
#' 
#' @param x a 'RasterBrick' or a three-dimensional array
#' @param depth depth map, generally a 'RasterLayer' object
#' @param layers vector of layer thickness
#' @param i0 a 'Raster' containing the number of soil laver just over the bedrock. Default is \code{NULL} and is then calculated.
#' @param verify logical. Default is \code{FALSE}. If it is \code{TRUE}, it verifies that function is working correctly. 
#' @param FUN function used for aggregation. If missing, \code{\link{identity}} is the default value.
#' @param divide.by.depth logical. If \code{TRUE} the function returns the 'mean' value, otherwise a a cumulate value. Default is \code{FALSE}. 
#' @param ... further argument for \code{FUN} 
#' 
#' @note \code{x} and \code{depth} or \code{i0} must cover the same spatial region.
#' 
#' @return a list of 'Raster' maps: 
#' 
#' \code{i0}   a 'Raster' containing the number of soil laver just over the bedrock
#' 
 
#' 
#' \code{z0}       a 'Raster' containing the depth of the center  of  the \code{i0}-th layer
#' 
#' \code{result}       a 'Raster' containing the aggregated map
#' 
#' @seealso  \code{\link{getvalues.brick.at.depth}},\code{\link{brick}}
#' @export
#' 
#' 
#' @examples
#' 
#' library(geotopbricks)
#' # The examples is the following R script conteined 
#' # in a 'inst' directory of the package source
#' f <- system.file("doc/examples/example.vertical.aggregate.brick.within.depth.R",
#' package="geotopbricks")
#' #  source(f) # Uncomment this line to run the example. 
#' # You can copy the example file using file.copy(from=f,to=....,...) See file.copy documentation
#' 
# library(rgdal)
# library(raster)
# library(zoo)
# library(methods)
# library(geotopbricks)
# library(soilwater)
#

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
# 
# ##  Water  column content: the soil water volume per unit area (along the hillslope), 'swc' function is used (see its documentation)  
# hw0 <- vertical.aggregate.brick.within.depth(psi,depth=bedrock_depth,layers=layers,verify=FALSE,FUN=swc,alpha=alphaVG,n=nVG,theta_sat=theta_sat,theta_res=theta_res) 
# ## expressed in millimeters
# # plot(hwo$result) # expressed in millimeters !!!! - Do Not Run!!!
# ##  Water column content: the soil water volume per unit topograhic area g the hillslope)
# hw <- hw0$result*cos(slope/180*pi)
# # plot(hw) # expressed in millimeters !!!! - Do Not Run!!!
#
#' 
#' 
#' 
#' 
#' 
#' 



vertical.aggregate.brick.within.depth<- function (x,depth=NULL,layers=NULL,i0=NULL,verify=FALSE,FUN=identity,divide.by.depth=FALSE,...) {
	

	
	if (is.null(layers)) layers <- array(1,nlayers(x))/nlayers(x)
		
	if (is.null(depth))  depth <- subset(x,1)*0+sum(layers)+1.0
	
	if (nlayers(x)!=length(layers)) {
		
		print("Error in vertical.aggregate.brick.within.depth: dimensions mismatching between x and layers ")
		return(-1)
	}
	
	if (nrow(x)!=nrow(depth)) {
		
		print("Error in vertical.aggregate.brick.within.depth: dimensions (rows) mismatching between x and depth ")
		return(-1)
	}
	
	if (ncol(x)!=ncol(depth)) {
		
		print("Error in vertical.aggregate.brick.within.depth: dimensions (columns) mismatching between x and depth ")
		return(-1)
	}
	

	
	z <- layers/2
	
	
	
	
	for (i in 2:length(z)) {
		
		z[i] <- z[i-1]+(layers[i-1]+layers[i])/2
	}
	
	
	if (is.null(i0)) {
		
		i0 <- depth*0+0

	
	
		for (i in 1:length(z)) {
			
			i0[depth>z[i]] <- i0[depth>z[i]]+1
		}



	} else if (nrow(x)!=nrow(i0)) {
		
		print("Error in vertical.aggregate.brick.within.depth: dimensions (rows) mismatching between x and i0 ")
		return(-1)
	} else if (ncol(x)!=ncol(i0)) {

		print("Error in vertical.aggregate.brick.within.depth: dimensions (columns) mismatching between x and i0 ")
		return(-1)
	} 
	
		
		
	i0[i0<1] <- 0 
	i0[i0>=length(z)] <- length(z)	
	zmax <- z[length(z)]+layers[length(z)]/2
	
	depth[depth>zmax] <- zmax
	
	
	out <- list()
	
	out$i0 <- i0
	out$z0 <- depth*0.0
	out$result <- depth*0.0 
	
	
	for (i in 1:length(layers)) {
	##	print(i)
		temp0 <- subset(x,i)
		val_temp <- FUN(temp0,...)
		out$result[out$i0>i] <-  out$result[out$i0>i]+val_temp[out$i0>i]*layers[i]
		out$result[out$i0==i] <-  out$result[out$i0==i]+val_temp[out$i0==i]*(depth[out$i0==i]-z[i]+layers[i]/2)
		out$z0[out$i0==i] <- z[i]	
	}	
	
	if (divide.by.depth) out$result <- out$result/depth
	
#	
#	for (r in 1:nrow(i0)) {
#		for (c in 1:ncol(i0)) {
#	
#			kindex <- i0[r,c]
#			if (!is.na(kindex)) {
#			
#			
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
			
			 print("Error in vertical.aggregate.brick.within.depth: depth bedrock is not included between z0 and z1, function does not work correctly!")
			 return(-1)
			 
		 }
	}
	return(out)
}


