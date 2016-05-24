#'Deriving the ground-projected canopy area of individual trees using LiDAR-derived derived Canopy Height Model (CHM) 
#'
#'@description Compute and export the ground-projected areas of individual tree canopies detected within the LiDAR-derived Canopy Height Model (CHM) 
#'
#'@usage ForestCAS(chm, loc, maxcrown, exclusion)
#'
#'@param chm A LiDAR-derived Canopy Height Model (CHM) RasterLayer or SpatialGridDataFrame file.
#'@param loc A matrix or dataframe with three columns (tree xy coordinates and height).
#'@param maxcrown A single value of the maximum individual tree crown radius expected. Default 10.0 m.
#'@param exclusion A single value from 0 to 1 that represents the \out{\%} of pixel exclusion. E.g. a value of 0.5 will exclude all of the pixels for a single tree that has a height value of less than 50\out{\%} of the maximum height from the same tree. Default value is 0.3. 
#'@return Returns a list that contains the individual tree canopy boundary polygons and the 4-column matrix with the tree xy coordinates, heights and ground-projected canopy area (with units of square meters).  
#'@author Carlos Alberto Silva
#'@examples
#'\dontrun{
#'
#'# Import the LiDAR-derived CHM file
#'data(chm) # or set a CHM. e.g. chm<-raster("CHM_stand.asc") 
#'
#'# Set the loc parameter
#'sCHM<-CHMsmoothing(chm, filter="mean", ws=5) # smoothing CHM
#'loc<-FindTreesCHM(sCHM, fws=5, minht=8)      # or import a tree list
#'
#'# Set the maxcrown parameter
#'maxcrown=10.0 
#'
#'# Set the exclusion parameter
#'exclusion=0.3 # 30
#'
#'# Compute individual tree detection canopy area
#'canopy<-ForestCAS(chm, loc, maxcrown, exclusion)
#'
#'#==================================================================================#
#'# Retrieving the boundary for individual tree detection and canopy area calculation
#'#==================================================================================#
#'boundaryTrees<-canopy[[1]]
#
#'# Plotting the individual tree canopy boundary over the CHM
#'plot(chm, main="LiDAR-derived CHM") 
#'plot(boundaryTrees, add=T, border='red', bg='transparent') # adding tree canopy boundary
#'
#'#============================================================================#
#'# Retrieving the list of individual trees detected for canopy area calculation
#'#============================================================================#
#'canopyList<-canopy[[2]] # list of ground-projected areas of individual tree canopies
#'summary(canopyList)     # summary 
#'
#'# Spatial location of the trees
#'library(sp)
#'XY<-SpatialPoints(canopyList[,1:2])    # Spatial points
#'plot(XY, col="black", add=T, pch="*")  # adding tree location to the plot
#'} 
#'@importFrom spatstat disc
#'@importFrom sp coordinates gridded Polygon Polygons SpatialPolygons over
#'@importFrom sp over SpatialPoints SpatialPolygonsDataFrame SpatialGridDataFrame SpatialPixelsDataFrame
#'@importFrom deldir tile.list deldir 
#'@importFrom plyr ddply
#'@importFrom raster raster rasterToPolygons boundaries
#'@export
ForestCAS<-function(chm,loc,maxcrown=10,exclusion=0.3) {

  if (class(chm)!="RasterLayer" & class(chm)!="SpatialGridDataFrame") {stop("The chm is invalid. It must to be a RasterLayer or SpatialGridDataFrame'")}
  if (ncol(loc)!=3) {stop("The input loc is invalid. It must to be 3-column matrix or dataframe with the x,y coordinates and heights of the individual trees")}
  if (class(maxcrown)!="numeric") {stop("The maxcrown parameter is invalid. It is not a numeric input")}
  if (class(exclusion)!="numeric") {stop("The exclusion parameter is invalid. It is not a numeric input")}
  if (exclusion >=1) {stop("The exclusion parameter is invalid. It must to be less than 1numeric input")}

  if (class(chm)=="RasterLayer"){ chm<-as(chm, "SpatialGridDataFrame")}
    
  Hthreshold<-min(loc[,3])*exclusion
  polys<-list() 
  width<-numeric()  
  
  DF2raster<-function(h.mH, i){
    
    h.mHkl<-subset(h.mH,h.mH[,4]==levels(factor(h.mH[,4]))[i])
    
    if (nrow(h.mHkl==1)==TRUE) {
      h.mHkl<-rbind(h.mHkl,c(h.mHkl[,1],h.mHkl[,2]+0.005,h.mHkl[,3]+0.005,h.mHkl[,4]))}
    
    spP <- cbind(h.mHkl[,2:3],h.mHkl[,1],h.mHkl[,4])
    colnames(spP)<-c("x","y","z","id")
    #coordinates(spP)<- c("x", "y")
    #suppressWarnings(gridded(spP) <- TRUE)
    m <- suppressWarnings(SpatialPixelsDataFrame(points=spP[c("x", "y")], data = spP))
    
    rasterDF <- raster(m)
    hhg<- boundaries(rasterDF, type='outer') 
    p <- rasterToPolygons(hhg, dissolve=TRUE)
    sp.polys <- p[1,]
    return(sp.polys)
  }
  
  for(i in 1:nrow(loc)) { 
    width[i]=maxcrown
    discbuff<-disc(radius=width[i], centre=c(loc$x[i], loc$y[i])) 
    discpoly<-Polygon(rbind(cbind(discbuff$bdry[[1]]$x, 
                                  y=discbuff$bdry[[1]]$y), c(discbuff$bdry[[1]]$x[1], 
                                                             y=discbuff$bdry[[1]]$y[1]))) 
    polys<-c(polys, discpoly) 
  } 
  
  spolys<-list() 
  for(i in 1:length(polys)) { 
    spolybuff<-Polygons(list(polys[[i]]), ID=row.names(loc)[i]) 
    spolys<-c(spolys, spolybuff) 
  } 
  polybuffs<-SpatialPolygons(spolys) 

  chmdf<-as.data.frame(chm,xy=TRUE)
  Points.Ply<-over(SpatialPoints(chmdf[,2:3]),polybuffs) 
  Points.PlyD<-cbind(chmdf,Points.Ply) 
  Points.PlyD<-na.omit(Points.PlyD) 
  vor =  deldir(loc[,1], loc[,2], z=loc[,3],suppressMsge=T)
  tile = tile.list(vor)
  polys = vector(mode='list', length=length(tile))
  
 for (i in seq(along=polys)) {
    pcrds = cbind(tile[[i]]$x, tile[[i]]$y)
    pcrds = rbind(pcrds, pcrds[1,])
    polys[[i]] = Polygons(list(Polygon(pcrds)), ID=as.character(i))
  }
  
  SP = SpatialPolygons(polys)
  veronoi = SpatialPolygonsDataFrame(SP, data=data.frame(x=loc[,1], 
                                                         y=loc[,2], row.names=sapply(slot(SP, 'polygons'), 
                                                                                          function(x) slot(x, 'ID'))))
  chmV<-over(SpatialPoints(Points.PlyD[,2:3]), SP)
  RpD<-cbind(Points.PlyD,chmV)
  RpD.filter<-subset(RpD[,1:5],RpD[,1]>=Hthreshold)
  RpD.filter<-cbind(RpD.filter[,1:3],RpD.filter[,5])
  colnames(RpD.filter)<-c("z","x","y","g")
  h.mH<-ddply(RpD.filter,"g", function (RpD.filter)
    subset(RpD.filter,RpD.filter[,1]>= max(RpD.filter[,1])*exclusion))
   
  for ( j in 1:nlevels(factor(h.mH[,4]))){
    assign(paste0("SP.polys", j), DF2raster(h.mH,j))
    cat (".");flush.console()}
  
  polygons <- slot(get("SP.polys1"), "polygons")
  
  for (i in 1:nlevels(factor(h.mH[,4]))) {
    data.loc <- get(paste0("SP.polys",i))
    polygons <- c(slot(data.loc, "polygons"),polygons)
  }
  
  for (i in 1:length(polygons)) {
    slot(polygons[[i]], "ID") <- paste(i)
  }
  
  spatialPolygons <- SpatialPolygons(polygons)
  spdf <- SpatialPolygonsDataFrame(spatialPolygons, 
                                   data.frame(Trees=1:length(polygons)))
  
  options(scipen=10)
  spdf<-spdf[spdf@data[-length(polygons),],]
  areaList<-sapply(slot(spdf, "polygons"), slot, "area")
  canopyTable<-cbind(loc,areaList)
  colnames(canopyTable)<-c("x","y","z","ca")
  result=list(spdf,canopyTable) 
  return(result)
  
}