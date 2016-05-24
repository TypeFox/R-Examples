#'3D convex hull of the individual tree LiDAR-derived point cloud 
#'
#'@description Compute and plot the 3D convex hull (and its surface area and volume) of the individual tree LiDAR-derived 
#'point cloud.
#'
#'@usage chullLiDAR3D(xyzid,plotit=TRUE,col="forestgreen",alpha=0.8)
#'
#'@param xyzid A matrix with four columns (xyz coordinates and tree id).
#'@param plotit Logical. If FALSE, returns only volume and surface area.
#'@param col A vector or a caracter of the convex hull color.
#'@param alpha A vector or a caracter of the convex hull transparency (0-1).
#'@return A list with components 'crownvolume' and 'crownsurface', giving the
#'volume and surface of the convex hull.
#'@author Carlos Alberto Silva. Uses code by Remko Duursma (\emph{YplantQMC} package,see "crownhull"). 
#'@references \url{www.qhull.org} and \emph{geometry} package (see \code{\link[geometry]{convhulln}}).
#'
#'@examples
#'
#'# Importing LAS file:
#'LASfile <- system.file("extdata", "LASexample1.las", package="rLiDAR")
#'
#'# Reading LAS file
#'LAS<-readLAS(LASfile,short=TRUE)
#'
#'# Setring the xyz coordinates and subsetting the data
#'xyz<-subset(LAS[,1:3],LAS[,3] >= 1.37)
#'
#'# Finding clusters
#'clLAS<-kmeans(xyz, 32)
#'
#'# Set the id vector
#'id<-as.factor(clLAS$cluster)
#'
#'#=================================================#
#'# Example 01
#'#=================================================#
#'# Set the alpha
#'alpha<-0.6
#'
#'# Set the plotCAS parameter
#'plotit=TRUE
#'
#'# Set the convex hull color
#'col="forestgreen"
#' 
#'# Combining xyz and id
#'xyzid<-cbind(xyz,id)
#' 
#'# Get the volume and surface area
#'library(rgl)
#'open3d() 
#'volumeList<-chullLiDAR3D(xyzid=xyzid, plotit=plotit, col=col,alpha=alpha)
#'summary(volumeList) # summary
#'
#'plot3d(xyzid[,1:3], add=TRUE)   # add the 3D point cloud
#'axes3d(c("x+", "y-", "z-"))                 # axes
#'grid3d(side=c('x+', 'y-', 'z'), col="gray") # grid
#'title3d(xlab = "UTM Easthing", ylab = "UTM Northing",zlab = "Height", col="red")
#'aspect3d(1,1,0.7) # scale
#'
#'#=================================================#
#'# Example 02
#'#=================================================#
#'# Set the alpha
#'alpha<-0.85
#'
#'# Set the plotCAS parameter
#'plotit=TRUE
#'
#'# Set the convex hull color
#'col=levels(factor(id))
#' 
#'# Combining xyz and id
#'xyzid<-cbind(xyz,id)
#' 
#'# Get the volume and surface area
#'open3d() 
#'volumeList<-chullLiDAR3D(xyzid=xyzid, plotit=plotit, col=col,alpha=alpha)
#'summary(volumeList)
#'
#'# Add other plot parameters
#'plot3d(xyzid[,1:3], col=xyzid[,4], add=TRUE)   # add the 3D point cloud
#'axes3d(c("x+", "y-", "z-"))                 # axes
#'grid3d(side=c('x+', 'y-', 'z'), col="gray") # grid
#'title3d(xlab = "UTM Easthing", ylab = "UTM Northing",zlab = "Height", col="red")
#'aspect3d(1,1,0.7) # scale
#'
#'@importFrom rgl triangles3d
#'@importFrom geometry convhulln
#'@export
chullLiDAR3D<-function(xyzid,plotit=TRUE,col="forestgreen", alpha=0.8) {

  VolumeList<-matrix(,ncol=3)[-1,]
  nlevels<-as.numeric(levels(factor(xyzid[,4])))
  
  for ( i in nlevels) {
    
    xyz<-subset(xyzid[,1:3],xyzid[,4]==i)
    
    if (nrow(xyz) <= 4) {crownvolume=0;crownsurface=0} else {
      
      
      if (plotit==TRUE) {
        if (length(col) == 1) {col2=col} 
          else { 
            if (length(nlevels)==length(col)) {col2=col[i]} 
              else { stop(cat("The col parameter doesn't has the same the lenght number of trees.")) }}
        if (length(alpha) == 1) {alpha2=alpha} 
            else { 
              if (length(nlevels)==length(alpha)) {alpha2=alpha[i]} 
                else { stop(cat("The alpha parameter doesn't has the same the lenght of number of trees.")) }}
              
        cat (".");flush.console()
        volume<-crownhull(xyz,plotit=TRUE,col=col2, alpha=alpha2)
      
      } else {
      volume<-crownhull(xyz,plotit=FALSE)        
      }  
      crownvolume<-as.numeric(volume[1])
      crownsurface<-as.numeric(volume[2]) }
    
    VolumeList<-data.frame(rbind(VolumeList,cbind(i,crownvolume,crownsurface)))
  }
  colnames(VolumeList)<-c("Tree", "crownvolume", "crownsurface")
  return(VolumeList)
}

