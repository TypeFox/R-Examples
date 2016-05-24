#'@name pathdistGen
#'@title Generate a stack of path distance raster objects
#'@description Generate a stack of path accumulated distance raster objects
#'@author Joseph Stachelek
#'@param spdf SpatialPointsDataFrame object
#'@param costras RasterLayer cost raster
#'@param range numeric. Range of interpolation neighborhood
#'@param yearmon character. String specifying the name of the spdf
#'@param progressbar logical show progressbar during processing?
#'@return RasterStack object of path distances
#'@import raster
#'@import gdistance
#'@importFrom utils setTxtProgressBar
#'@export
#'@examples
#'spdf <- data.frame(rnorm(2))
#'xy <- data.frame(x = c(4, 2), y = c(8, 4))
#'coordinates(spdf) <- xy
#'
#'m <- matrix(NA, 10, 10)
#'costras <- raster(m, xmn = 0, xmx = ncol(m), ymn = 0, ymx = nrow(m))
#'costras[] <- runif(ncell(costras), min = 1, max = 10)
#'#introduce spatial gradient
#'for(i in 1:nrow(costras)){
#'  costras[i,] <- costras[i,] + i
#'  costras[,i] <- costras[,i] + i
#'}
#'
#'rstack <- pathdistGen(spdf, costras, 100)

'pathdistGen' <- function(spdf, costras, range, yearmon = "default", progressbar = TRUE){
  
  if(class(spdf) != "SpatialPointsDataFrame"){
    stop("spdf object must be of class SpatialPointsDataFrame")
  }
  
  ipdw_range <- range / raster::res(costras)[1] / 2 #this is a per cell distance
  
  #start interpolation
  #calculate conductances hence 1/max(x)
  trans <- gdistance::transition(costras, function(x) 1 / max(x), directions = 16)
  i <- 1
  coord <- spdf[i,]
  costsurf <- gdistance::accCost(trans, coord)
  
  ipdw_dist <- hist(costsurf, plot = F)$breaks[2]
  if(ipdw_dist < ipdw_range){
    ipdw_dist <- ipdw_range    
  }
  
  if(progressbar == TRUE){pb <- utils::txtProgressBar(max = nrow(spdf), style = 3)}
  
    for(i in 1:nrow(spdf)){
      coord <- spdf[i,]
      costsurf <- gdistance::accCost(trans, coord)
      costsurf_reclass <- raster::reclassify(costsurf, c(ipdw_dist, +Inf, NA, ipdw_range, ipdw_dist, ipdw_range)) #raster cells are 1 map unit
      costsurf_scaled <- ((ipdw_range/costsurf_reclass)^2)
      costsurf_scaled_reclass <- raster::reclassify(costsurf_scaled, c(-Inf, 1, 0))
      
      #check output with - zoom(A4,breaks=seq(from=0,to=5,by=1),col=rainbow(5))
      #showTmpFiles()
      
      raster::writeRaster(costsurf_scaled_reclass, filename = file.path(tempdir(), paste(yearmon, "A4ras", i, ".grd", sep = "")), overwrite = T, NAflag = -99999)
      if(progressbar == TRUE){setTxtProgressBar(pb, i)}      
    }
  if(progressbar == TRUE){close(pb)}
     
  
  #create raster stack
  raster_flist <- list.files(path = file.path(tempdir()), pattern = paste(yearmon, "A4ras*", sep = ""), full.names = T)
  raster_flist <- raster_flist[grep(".grd", raster_flist, fixed = T)]
  as.numeric(gsub('.*A4ras([0123456789]*)\\.grd$', '\\1', raster_flist)) -> fileNum
  raster_flist <- raster_flist[order(fileNum)]
  rstack <- raster::stack(raster_flist)
  rstack <- raster::reclassify(rstack, cbind(-99999, NA))
  file.remove(list.files(path = file.path(tempdir()), pattern = paste(yearmon, "A4ras*", sep = ""), full.names = T))
  
  return(rstack)
}
