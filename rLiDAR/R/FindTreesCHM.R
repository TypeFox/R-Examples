#'Individual tree detection whitin the LiDAR-derived Canopy Height Model (CHM) 
#'
#'@description Detects and computes the location and height of individual trees within the LiDAR-derived Canopy Height Model (CHM). The algorithm implemented in this function is local maximum with a fixed window size. 
#'
#'@usage FindTreesCHM(chm,fws,minht)
#'
#'@param chm A LiDAR-derived Canopy Height Model (CHM) raster  file.
#'@param fws A single dimension (in raster grid cell units) of fixed square window size, e.g. 3, 5, 7 and so on. Default is 5. 
#'@param minht Height threshold. Detect individual trees above specified height threshold, e.g. 1.37, 2.0, 3.5 m and so on. Default is 1.37 m.
#'@return Returns A matrix with four columns (tree id, xy coordinates, and height).
#'@author Carlos Alberto Silva
#'@examples
#'
#'# Importing the LiDAR-derived CHM raster file
#'data(chm) # or set a CHM. e.g. chm<-raster("CHM_stand.asc") 
#'
#'# Smoothing CHM
#'schm<-CHMsmoothing(chm, "mean", 5)
#'
#'# Setting the fws:
#'fws<-5 # dimention 5x5
#'
#'# Setting the specified height above ground for detectionbreak
#'minht<-8.0
#'
#'# Getting the individual tree detection list
#'treeList<-FindTreesCHM(schm, fws, minht)
#'summary(treeList)
#'
#'# Plotting the individual tree location on the CHM
#'library(raster)
#'plot(chm, main="LiDAR-derived CHM")
#'library(sp) 
#'XY<-SpatialPoints(treeList[,1:2]) # Spatial points
#'plot(XY, add=TRUE, col="red")        # plotthing tree location
#'
#'@importFrom raster raster focal xyFromCell Which projection
#'@importFrom sp SpatialPoints over SpatialGridDataFrame 
#'@export
FindTreesCHM<-function(chm, fws=5,minht=1.37) {
  
  if (class(chm)[1]!='RasterLayer') {chm<-raster(chm)}
  if (class(fws)!="numeric") {stop("The fws parameter is invalid. It is not a numeric input")}
  if (class(minht)!="numeric") {stop("The minht parameter is invalid. It is not a numeric input")}
  
  w<-matrix(c(rep(1,fws*fws)),nrow=fws,ncol=fws)
    
  chm[chm < minht]<-NA
  
  f <- function(chm) max(chm)
  
  rlocalmax <- focal(chm, fun=f, w=w, pad=TRUE, padValue=NA)
  
  setNull<- chm==rlocalmax
  XYmax <- SpatialPoints(xyFromCell(setNull, Which(setNull==1, cells=TRUE)))
  
  #proj<-projection(chm)
  #projection(XYmax)=proj
  htExtract<-over(XYmax,as(chm, "SpatialGridDataFrame"))
  treeList<-cbind(XYmax,htExtract)
  
  colnames(treeList)<-c("x","y","height")
    
  return(treeList)
 
}
