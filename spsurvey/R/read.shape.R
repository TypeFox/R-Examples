read.shape <- function(filename=NULL) {

################################################################################
# Function: read.shape
# Purpose: Read an ESRI shapefile
# Programmer: Tom Kincaid
# Date: March 1, 2005
# Last Revised: October 8, 2008
# Description:
#   This function reads either a single shapefile or multiple shapefiles.  For 
#   multiple shapefiles, all of the shapefiles must be the same type, i.e., 
#   point, polyline, or polygon.
# Arguments:
#   filename = name of the shapefile without any extension.  If filename equals
#     a shapefile name, then that shapefile is read.  If filename equals NULL,
#     then all of the shapefiles in the working directory are read.  The default
#     is NULL.
# Results:
#   An sp package object containing information in the shapefile.  The object is
#   assigned class "SpatialPointsDataFrame", "SpatialLinesDataFrame", or
#   "SpatialPolygonsDataFrame" corresponding to the shapefile type, i.e., point,
#   polyline, or polygon, respectively.  For further information regarding the
#   output object, see documentation for the sp package.
# Other Functions Required:
#   readShapeFile - C function to read a single shapefile or multiple shapefiles
#   SpatialPoints - sp package function to create an object of class
#     SpatialPoints
#   SpatialPointsDataFrame - sp package function to create an object of class
#     SpatialPointsDataFrame
#   shape2spList - function to create an object of class Lines for a lines
#      shapefile or class Polygons for a polygons shapefile
#   SpatialLines - sp package function to create an object of class SpatialLines
#   SpatialLinesDataFrame - sp package function to create an object of class
#     SpatialLinesDataFrame
#   SpatialPolygons - sp package function to create an object of class
#     SpatialPolygons
#   SpatialPolygonsDataFrame - sp package function to create an object of class
#     SpatialPolygonsDataFrame
################################################################################

# If necessary, strip the file extension from the file name

   if(!is.null(filename)) {
      nc <- nchar(filename)
      if(substr(filename, nc-3, nc) == ".shp") {
         filename <- substr(filename, 1, nc-4)
      }
   }

# Read the shapefile

   sfile <- .Call("readShapeFile", filename)
   if(is.null(sfile[[1]]))
      stop("\nAn error occurred while reading the shapefile(s) in the working directory.")

# Convert character vectors to factors in the attributes data frame

   ind <- sapply(sfile$att.data, is.character)
   if(any(ind)) {
      for(i in (1:length(sfile$att.data))[ind])
         sfile$att.data[,i] <- as.factor(sfile$att.data[,i])
   }

# Create an sp package object

   att.data <- sfile$att.data
   shapes <- sfile$Shapes
   n <- length(shapes)
   IDs <- as.character(1:n)
   rownames(att.data) <- IDs
   shp.type <- attr(sfile$Shapes, "shp.type")
   if(shp.type == "point") {
      SpointsMat <- matrix(0, nrow=n, ncol=2)
      rownames(SpointsMat) <- IDs
      for(i in 1:n) {
        SpointsMat[i,] <- shapes[[i]]$verts
      }
      sp.obj <- SpatialPointsDataFrame(coords=SpatialPoints(coords=SpointsMat),
         data=att.data)
   } else if(shp.type == "arc") {
      SlinesList <- vector(mode="list", length=n)
      for(i in 1:n) {
        SlinesList[[i]] <- shape2spList(shapes[[i]], shp.type, IDs[i])
      }
      sp.obj <- SpatialLinesDataFrame(sl=SpatialLines(LinesList=SlinesList),
         data=att.data)
#      for(i in 1:n) {
#         sp.obj@lines[[i]]$length <- shapes[[i]]$length
#      }
   } else if(shp.type == "poly"){
      PolygonsList <- vector(mode="list", length=n)
      for(i in 1:n) {
        PolygonsList[[i]] <- shape2spList(shapes[[i]], shp.type, IDs[i])
      }
      sp.obj <- SpatialPolygonsDataFrame(Sr=SpatialPolygons(Srl=PolygonsList),
         data=att.data)
      for(i in 1:n) {
         sp.obj@polygons[[i]]@area <- sum(shapes[[i]]$areas)
         nParts <- length(sp.obj@polygons[[i]]@Polygons)
         for(j in 1:nParts) {
            sp.obj@polygons[[i]]@Polygons[[j]]@area <- shapes[[i]]$areas[j]
         }
      }
   } else {
      stop(paste("\nShapefile type", shp.type, "is not recognized."))
   }

# Return the sp package object

   sp.obj
}
   