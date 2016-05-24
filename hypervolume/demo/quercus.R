if (exists('doHypervolumeQuercusDemo')==TRUE)
{
  require(raster)
  require(maps)
  
  # load in lat/lon data
  data('quercus') 
  data_alba = subset(quercus, Species=="Quercus alba")[,c("Longitude","Latitude")]
  data_rubra = subset(quercus, Species=="Quercus rubra")[,c("Longitude","Latitude")]
  
  # get worldclim data from internet
  climatelayers = getData('worldclim', var='bio', res=10, path=tempdir())
  
  # z-transform climate layers to make axes comparable
  climatelayers_ss = climatelayers[[c(1,4,12,15)]]
  for (i in 1:nlayers(climatelayers_ss))
  {
    climatelayers_ss[[i]] <- (climatelayers_ss[[i]] - cellStats(climatelayers_ss[[i]], 'mean')) / cellStats(climatelayers_ss[[i]], 'sd') 
  }
  
  # extract transformed climate values
  climate_alba = extract(climatelayers_ss, data_alba)
  climate_rubra = extract(climatelayers_ss, data_rubra)
  
  # compute hypervolumes with auto-bandwidth for both species
  hv_alba = hypervolume(climate_alba,quantile=0.0,reps=1000,bandwidth=estimate_bandwidth(climate_alba),name='alba')
  hv_rubra = hypervolume(climate_rubra,quantile=0.0,reps=1000,bandwidth=estimate_bandwidth(climate_rubra),name='rubra')
  
  # determine intersection and unique components of the overlap
  hv_set = hypervolume_set(hv_alba, hv_rubra, check_memory=FALSE)
  
  # put all the output volumes in one convenient place
  volumes <- c(Alba=get_volume(hv_alba), Rubra=get_volume(hv_rubra), get_volume(hv_set))
  
  # do species distribution modeling
  # get all the climate values
  climatevalues = data.frame(getValues(climatelayers_ss))
  
  rubra_inout = hypervolume_inclusion_test(hv_rubra, climatevalues)
  alba_inout = hypervolume_inclusion_test(hv_alba, climatevalues)
  
  # convert to rasters by setting values in rasters with same extent/resolution
  rubra_map = raster(climatelayers_ss[[1]]); values(rubra_map) <- rubra_inout
  alba_map = raster(climatelayers_ss[[1]]); values(alba_map) <- alba_inout
  
  # then barplot of hypervolumes of each component
  barplot(volumes,horiz=TRUE,las=2,main="Hypervolume")
  
  # then pairs plot of the set operations
  plot(hv_set)
  
  # plot the geographic projections of the ranges
  par(mfrow=c(1,2))
  plot(rubra_map,col=c(NA,rgb(1,0,0,0.5)),legend=FALSE,xlim=c(-100,-50),ylim=c(20,60),main='Quercus rubra')
  map('world',add=TRUE)
  points(Latitude~Longitude,data=data_rubra,pch=3,cex=0.25)
  
  plot(alba_map,col=c(NA,rgb(0,0,1,0.5)),legend=FALSE,xlim=c(-100,-50),ylim=c(20,60),main='Quercus alba')
  map('world',add=TRUE)
  points(Latitude~Longitude,data=data_alba,pch=3,cex=0.25)
  par(mfrow=c(1,1))
  
  rm(doHypervolumeQuercusDemo)
} else
{
  message('Demo does not run by default to meet CRAN runtime requirements.')
  message('This demo requires internet access to download 10MB of climate data and')
  message('will take approximately 3 minutes to run.')
  message('To run the demo, type')
  message('\tdoHypervolumeQuercusDemo=TRUE')
  message('\tdemo(quercus)')
  message('at the R command line prompt.')
}