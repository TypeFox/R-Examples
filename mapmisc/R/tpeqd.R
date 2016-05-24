tpeqdProj4string = function(
    lon1, lat1, lon2, lat2, 
    x=0,y=0, datum='WGS84',
    ellps='WGS84', units='m',
    crs=TRUE) {
  
  res = paste(
      "+proj=tpeqd",
      " +lat_1=",lat1,
      " +lon_1=",lon1,
      " +lat_2=",lat2,
      " +lon_2=",lon2,
      " +x_0=",x,
      " +y_0=",y,
      " +ellps=",ellps,
      " +datum=",datum,
      " +units=",units,
      " +no_defs",sep=''
      )
  if(crs){
    res = lapply(res, CRS)
  }
  
}

tpeqd = function(x, offset=c(0,0)){
  
  if(length(grep("^SpatialPoints", class(x)))){
    if(!isLonLat(crs(x)) & 
        requireNamespace('rgdal', quietly=TRUE)) {
      x = spTransform(x, crsLL)
    }
    x = coordinates(x)
  }
  
  x = as.matrix(x)[1:2,1:2]
  x = x[order(x[,2],decreasing=TRUE),]
  
  
  res= tpeqdProj4string(
      x[1,1],x[1,2],x[2,1],x[2,2],
      x=offset[1],y=offset[2]
      )
  if(length(res)[[1]])
    res = res[[1]]
      
  res
  
}