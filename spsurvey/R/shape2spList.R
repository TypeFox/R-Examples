shape2spList <- function (shape, shp.type, ID) {

################################################################################
# Function: shape2spList
# Purpose: Create an object of class Lines or class Polygons
# Programmer: Tom Kincaid
# Date: July 21, 2005
# Last Revised: January 18, 2006
# Description:
#   This function creates an object of class Lines for a Polyline shapefile or
#   class Polygons for a Polygon shapefile.
# Arguments:
#   shape = a single record from the .shp file of the shapefile.
#   shp.type - the type of shapefile, which is either "arc" for a Polyline
#      shapefile or "poly" for a Polygon shapefile.
#   ID - the shape ID value, i.e., the shapefile record number.
# Results:
#   An object of class Lines for a Polyline shapefile or class Polygons for a
#   Polygon shapefile - see documentation for the sp package for further
#   details.
# Other Functions Required:
#   Line - sp package function to create an object of class Line
#   Lines - sp package function to create an object of class Lines
#   Polygon - sp package function to create an object of class Polygon
#   Polygons - sp package function to create an object of class Polygons
################################################################################

   nParts <- shape$nParts
   nVerts <- shape$nVerts
   Pstart <- shape$Pstart
   from <- integer(nParts)
   to <- integer(nParts)
   from[1] <- 1
   for(j in 1:nParts) {
      if(j == nParts) {
         to[j] <- nVerts
      } else {
         to[j] <- Pstart[j + 1]
         from[j + 1] <- to[j] + 1
      }
   }
   temp <- vector(mode="list", length=nParts)
   if(shp.type == "arc") {
      for(i in 1:nParts) {
         temp[[i]] <- Line(coords=shape$verts[from[i]:to[i],])
      }
      Lines <- Lines(slinelist=temp, ID=ID)
      return(Lines)
   } else {
      hole.ind <- as.logical(attr(shape, "RingDir") - 1)
      for(i in 1:nParts) {
         temp[[i]] <- Polygon(coords=shape$verts[from[i]:to[i],],
            hole=hole.ind[i])
      }
      Polygons <- Polygons(srl=temp, ID=ID)
      return(Polygons)
   }
}
