#'2D Convex hull of individual tree LiDAR-derived point cloud 
#'
#'@description Compute and plot the 2D convex hull of individual tree LiDAR-derived point cloud 
#'
#'@usage chullLiDAR2D(xyid)
#'
#'@param xyid A 3-column matrix with the x, y coordinates and points id of the LiDAR point cloud.
#'@return Returns A list with components "chullPolygon" and "chullArea", giving the polygon and area  of the convex hull.  
#'@author Carlos Alberto Silva
#'@references \emph{grDevices} package,see \code{\link[grDevices]{chull}}.
#'@examples
#'
#'# Importing LAS file:
#'LASfile <- system.file("extdata", "LASexample1.las", package="rLiDAR")
#'
#'# Reading LAS file
#'LAS<-readLAS(LASfile,short=TRUE)
#'
#'# Height subsetting the data
#'xyz<-subset(LAS[,1:3],LAS[,3] >= 1.37)
#'
#'# Getting LiDAR clusters
#'set.seed(1)
#'clLAS<-kmeans(xyz, 32)
#'
#'# Set the points id 
#'id<-as.factor(clLAS$cluster)
#'
#'# Set the xyid input
#'xyid<-cbind(xyz[,1:2],id)
#'
#'# Compute the LiDAR convex hull of the clusters 
#'chullTrees<-chullLiDAR2D(xyid)
#'
#'# Plotting the LiDAR convex hull
#'library(sp)
#'plot(SpatialPoints(xyid[,1:2]),cex=0.5,col=xyid[,3])
#'plot(chullTrees$chullPolygon,add=TRUE, border='green')
#'
#'# Get the ground-projected area of LiDAR convex hull
#'chullList<-chullTrees$chullArea 
#'summary(chullList) # summary 
#'
#'@export
chullLiDAR2D<-function(xyid) {
  
  spdfList<-list()
  for ( i in levels(factor(xyid[,3]))) {
    
    subSet<-subset(xyid,xyid[,3]==i)
    dat<-subSet[,1:2]
    ch <- chull(dat)
    coords <- dat[c(ch, ch[1]), ]  
    sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=i)))
    spdfList[i]<-SpatialPolygonsDataFrame(sp_poly, data=data.frame(subSet[1,1], row.names=row.names(sp_poly)))
    cat (".");flush.console()
  } 
    
  polygons <- slot(spdfList[[1]], "polygons")
  
  for (i in levels(factor(xyid[,3]))) {
    data.loc <- spdfList[[i]]
    polygons <- c(slot(data.loc, "polygons"),polygons)
  }
  
  for (i in 1:length(polygons)) {
    slot(polygons[[i]], "ID") <- paste(i)
    
  }
  
  spatialPolygons <- SpatialPolygons(polygons)
  spdf <- SpatialPolygonsDataFrame(spatialPolygons, 
                                   data.frame(TreeID=1:length(polygons)))
  
  options(scipen=4)
  spdf<-spdf[spdf@data[-length(polygons),],]
  areaList<-as.numeric(sapply(slot(spdf, "polygons"), slot, "area"))
  canopyTable<-data.frame(cbind(as.numeric(levels(factor(xyid[,3]))),rev(areaList)))
  spdf@data$TreeID<-canopyTable[,1]
  spdf@data$GPA<-canopyTable[,2]
  colnames(canopyTable)<-c("TreeID","GPA")
  return(list(chullPolygon=spdf,chullArea=canopyTable)) 
      
}

