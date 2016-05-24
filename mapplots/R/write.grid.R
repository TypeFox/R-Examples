write.grid <- function (grd, file, type='csv') 
{
  lon0 <- as.numeric(rownames(grd))
  lat0 <- as.numeric(colnames(grd))
  lon = rep(lon0, times = dim(grd)[2])
  lat = rep(lat0, each = dim(grd)[1])
  value = as.numeric(grd)
  df <- na.omit(data.frame(lon, lat, value)) 
  if(type=='shape'){
    if (requireNamespace("shapefiles", quietly = TRUE)==FALSE) {
      stop('The package shapefiles is not loaded')
    } else {
    byx <- diff(range(lon0)) / (dim(grd)[1] - 1)
    byy <- diff(range(lat0)) / (dim(grd)[2] - 1)
    n <- nrow(df)    
    shpTable <- data.frame(Id=1:n,X=df$lon-byx/2,Y=df$lat-byy/2)
    shpTable <- rbind(shpTable,data.frame(Id=1:n,X=df$lon+byx/2,Y=df$lat-byy/2))
    shpTable <- rbind(shpTable,data.frame(Id=1:n,X=df$lon+byx/2,Y=df$lat+byy/2))
    shpTable <- rbind(shpTable,data.frame(Id=1:n,X=df$lon-byx/2,Y=df$lat+byy/2))
    shpTable <- rbind(shpTable,data.frame(Id=1:n,X=df$lon-byx/2,Y=df$lat-byy/2))
    shpTable <- shpTable[order(shpTable$Id),]
    attTable <- data.frame(Id=1:n,value=df$value)
    shp <- shapefiles::convert.to.shapefile(shpTable,attTable,'Id',5)
    shapefiles::write.shapefile(shp,file,arcgis=T)
    proj <- 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]'
    write.table(proj,paste(file,'.prj',sep=''),quote=F,row.names=F,col.names=F)
  }} else write.csv(df, file, row.names = F)
  
}