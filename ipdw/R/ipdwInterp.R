#'@name ipdwInterp
#'@title Inverse Distance Weighting with custom distances
#'@description This function takes a rasterstack of pathdistances and generates surfaces by weighting parameter values by these distances
#'@author Joseph Stachelek
#'@param spdf SpatialPointsDataFrame object
#'@param rstack RasterStack of path distances
#'@param paramlist chacter. String representing parameter names
#'@param overlapped logical. Default is FALSE, specify TRUE if some points lie on top of barriers
#'@param yearmon character. String specifying the name of the spdf
#'@param removefile logical. Remove files after processing?
#'@return RasterLayer
#'@import raster
#'@import gdistance
#'@export
#'@examples
#'spdf <- data.frame(rnorm(2))
#'xy <- data.frame(x = c(4, 2), y = c(8, 4))
#'coordinates(spdf) <- xy
#'m <- matrix(NA, 10, 10)
#'costras <- raster(m, xmn = 0, xmx = ncol(m), ymn = 0, ymx = nrow(m))
#'#introduce spatial gradient
#'costras[] <- runif(ncell(costras), min = 1, max = 10)
#'for(i in 1:nrow(costras)){
#'  costras[i,] <- costras[i,] + i
#'  costras[,i] <- costras[,i] + i
#'}
#'
#'rstack <- pathdistGen(spdf, costras, 100)
#'final.raster <- ipdwInterp(spdf, rstack, paramlist = c("rnorm.2."), overlapped = TRUE)
#'plot(final.raster)

'ipdwInterp' <- function(spdf, rstack, paramlist, overlapped = FALSE, yearmon = "default", removefile = TRUE){
  
  for(k in 1:length(paramlist)){
  	points_layers <- rm_na_pointslayers(param_name = paramlist[k], spdf = spdf, rstack = rstack)
  	spdf <- points_layers$spdf
  	rstack <- points_layers$rstack
  	
    #if(overlapped == TRUE){ #need to set na.rm = TRUE if points are on land?
    rstack.sum <- raster::calc(rstack, fun = function(x){sum(x, na.rm = TRUE)})
    rstack.sum <- raster::reclassify(rstack.sum, cbind(0, NA))
      
    #calculate the weight of the individual rasters 
    for(i in 1:dim(rstack)[3]){
      ras.weight <- rstack[[i]] / rstack.sum
      param.value <- data.frame(spdf[i, paramlist[k]])
      param.value2 <- as.vector(unlist(param.value[1]))
      ras.mult <- ras.weight * param.value2
      
      rf <- raster::writeRaster(ras.mult, filename = file.path(tempdir(), paste(paramlist[k], "A5ras", i, ".grd", sep = "")), overwrite = T)
    }
    
    raster_data_full <- list.files(path = file.path(tempdir()), pattern = paste(paramlist[k], "A5ras*", sep = ""), full.names = T)
    raster_data <- raster_data_full[grep(".grd", raster_data_full, fixed = T)]
    as.numeric(gsub('.*A5ras([0123456789]*)\\.grd$', '\\1', raster_data)) -> fileNum
    raster_data <- raster_data[order(fileNum)]
    
    #sum rasters to get final surface
    rstack.mult <- raster::stack(raster_data)
    finalraster <- raster::calc(rstack.mult, fun = function(x){sum(x, na.rm = TRUE)})
    
    if(overlapped == TRUE){
    	finalraster <- raster::reclassify(finalraster, cbind(0, NA))
    }
     
    r <- raster::rasterize(spdf, rstack[[1]], paramlist[k])
    finalraster <- raster::cover(r, finalraster)
    
    file.remove(raster_data_full)
    return(finalraster)
  }
#optional removal of path distances

  if(removefile == TRUE){
    file.remove(list.files(path = file.path(tempdir()), pattern = paste(yearmon, "A4ras*", sep = "")))
    file.remove(list.files(path = file.path(tempdir()), pattern = paste(paramlist[k], "A5ras*", sep = ""), full.names = T))
  }
}

#'@name rm_na_pointslayers
#'@title Remove NA SpatialPointsDataFrame features and drop correspoding raster stack layers
#'@description Remove NA SpatialPointsDataFrame features and drop correspoding raster stack layers
#'@param param_name character name of data column
#'@param spdf SpatialPointsDataFrame object
#'@param rstack RasterStack or RasterBrick
#'@export

rm_na_pointslayers <- function(param_name, spdf, rstack){
	param_index_x <- which(names(spdf) == param_name)
	param_na_y <- which(is.na(spdf@data[,param_index_x]))
	
	if(length(param_na_y) > 0){
		spdf <- spdf[-which(is.na(spdf@data[,param_index_x])),]  
		rstack <- raster::dropLayer(rstack, param_na_y)
	}
	list(spdf = spdf, rstack = rstack)
}