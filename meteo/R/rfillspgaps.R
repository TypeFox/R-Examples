rfillspgaps <- function(rasterLayer,
                      maskPol=NULL,
                      nmax =50,
                      zcol=1,
                      ...
                      ){

#   if(class(rasterLayer)=="RasterLayer") { r=rasterLayer} else{r= raster(rasterLayer) }
  r=as(rasterLayer,'SpatialGridDataFrame')
  crs= r@proj4string
  
  if(!is.null(maskPol)){
    rr <- over( r, as(maskPol,"SpatialPolygons"),  returnList = FALSE)
    r[[zcol]][is.na(rr)] <- c(-5555)
    r[[zcol]][is.na(r[[zcol]])] <-c(-9999)
    r= rasterToPoints(raster(r[zcol]), spatial = T , fun= function(x){ x != -5555})
    rr=NULL
  }
    
r <- as( r, 'SpatialPointsDataFrame' )
r[is.na(r@data[,zcol]),zcol] <- c(-9999)


data <- r[(r@data[,zcol]!=c(-9999)),]
data@proj4string <- CRS(as.character(NA))
newdata=r[r@data[,zcol]==c(-9999),]
newdata@proj4string <- CRS(as.character(NA))

if(is.numeric(zcol)) {zcol= names(data)[zcol] }

 gaps = idw(as.formula(paste(zcol,"~1",sep="")), data,newdata ,nmax =nmax,...)
 gaps=gaps[,1]
#  names(gaps)= zcol
 r[(r@data[,zcol]==c(-9999)),] <- gaps@data[,1]

 r= rasterFromXYZ(as.data.frame(r[,zcol])[,c(2,3,1)] , crs=crs)
 
 if(class(rasterLayer)!="RasterLayer"){
   r= as(r,class(rasterLayer))
 }
  return(r)
}
  

