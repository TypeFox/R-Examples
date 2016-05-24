sbcframe <- function(shapefilename = NULL, spframe = NULL, nrows = 5,
   dxdy = TRUE) {

################################################################################
# Function: sbcframe
# Purpose: Calculate spatial balance grid cell extent and proportions for a
#          sample frame
# Programmer: Tom Kincaid
# Date: September 29, 2011
# Revised: April 17, 2015
# Description:      
#   This function calculates spatial balance grid cell extent and proportions
#   for the sample frame.  
# Arguments:
#   shapefilename = name of the input shapefile.  If shapefilename equals NULL,
#     then the shapefile or shapefiles in the working directory are used.  The
#     default is NULL.
#   spframe = an sp package object of class SpatialPointsDataFrame,
#     SpatialLinesDataFrame, or SpatialPolygonsDataFrame that contains the
#     survey design frame.  The default is NULL.
#   nrows = number of rows (and columns) for the grid of cells.  The default is
#     5.
#   dxdy = indicator for equal x-coordinate and y-coordinate grid cell
#     increments, where TRUE means the increments are equal and FALSE means the
#     increments are not equal.  The default is TRUE.
# Results: 
#   A list containing the following components: (1) extent - the frame extent
#   for each grid cell, (2) prop - the frame proportion for each grid cell,
#  (3) xmin - the grid x-coordinate minimum value, (4) xmax - the grid
#   x-coordinate maximum value, (5) ymin - the grid y-coordinate minimum value,
#  (6) ymax - the grid y-coordinate maximum value, (7) dx - the grid cell
#   x-coordinate increment value, (8) dy - the grid cell y-coordinate increment
#   value, (9) xc - the vector of grid cell x-coordinates, and (10) yc - the
#   vector of grid cell y-coordinates.
# Other Functions Required:
#   readShapeFile - C function to read a single shapefile or multiple shapefiles
#   readShapeFilePts - C function to read the shp file of a point shapefile and
#     return a data frame containing the x-coordinates and y-coordinates for
#     elements in the frame
#   cell.wt - calculates number of points in a cell for a points object
#   insideLinearGridCell - C function to determine ID value and clipped polyline
#     length for shapefile records contained in the selected grid cells
#   insideAreaGridCell - C function to determine ID value and clipped polygon
#     area for shapefile records contained in the selected grid cells
################################################################################

# Check that either a shapefile name of a survey design frame object was provided
   if(is.null(shapefilename) & is.null(spframe))
      stop("\nEither a shapefile name or a survey design frame object must be provided.")

# If necessary, strip the file extension from the shapefile name
   if(!is.null(shapefilename)) {
      nc <- nchar(shapefilename)
      if(substr(shapefilename, nc-3, nc) == ".shp") {
         shapefilename <- substr(shapefilename, 1, nc-4)
      }
   }
# If a survey design frame object was provided, then create a temporary
# shapefile
   if(!is.null(spframe)) {
      shapefilename <- "shapefile0202"
      sp2shape(spframe, shapefilename)
   }

# Read the shapefile
   sfile <- .Call("readShapeFile", shapefilename)
   if(is.null(sfile[[1]]))
      stop("\nAn error occurred while reading the shapefile(s) in the working directory.")

# Determine the type of shapefile
   shp.type <- attr(sfile$Shapes, "shp.type")

# Determine the number of records in the shapefile
   nshps <- attr(sfile$Shapes, "nshps")

# Calculate the x-coordinate and y-coordinate increment values and create the
# vectors of grid x-coordinates and y-coordinates
   minbb <- attr(sfile$Shapes, "minbb")
   maxbb <- attr(sfile$Shapes, "maxbb")
   xmin <- minbb[1]
   ymin <- minbb[2]
   xmax <- maxbb[1]
   ymax <- maxbb[2]
   if(dxdy) {
      gridExtent = max((xmax - xmin), (ymax - ymin))
      xmin = xmin - gridExtent * 0.001
      xmax = xmin + gridExtent * 1.002
      ymin = ymin - gridExtent * 0.001
      ymax = ymin + gridExtent * 1.002
   } else {
      gridExtent = xmax - xmin;
      xmin = xmin - gridExtent * 0.001
      xmax = xmin + gridExtent * 1.002
      gridExtent = ymax - ymin
      ymin = ymin - gridExtent * 0.001
      ymax = ymin + gridExtent * 1.002
   }
   dx <- (xmax - xmin)/nrows
   dy <- (ymax - ymin)/nrows
   xc <- seq(xmin, xmax, length=(nrows+1))[-1]
   xc <- rep(xc, nrows)
   yc <- seq(ymin, ymax, length=(nrows+1))[-1]
   yc <- rep(yc, rep(nrows, nrows))
   ncells <- length(xc)

# Calculate grid cell extent and proportion for a point shapefile
   if(shp.type == "point") {
      temp <- .Call("readShapeFilePts", shapefilename)
      ptsframe <- data.frame(x=temp$x, y=temp$y, mdm=1)
      extent <- sapply(1:ncells, cell.wt, xc, yc, dx, dy, ptsframe)
      prop <- extent/sum(extent)

# Calculate grid cell extent and proportion for a polyline shapefile
   } else if(shp.type == "arc") {
      extent <- numeric(ncells)
      temp <- .Call("insideLinearGridCell", shapefilename, 1:nshps, 1:ncells,
         xc, yc, dx, dy)
      temp <- tapply(temp$recordLength, temp$cellID, sum)
      extent[as.numeric(names(temp))] <- temp
      prop <- extent/sum(extent)

# Calculate grid cell extent and proportion for a polygon shapefile
   } else if(shp.type == "poly") {
      extent <- numeric(ncells)
      temp <- .Call("insideAreaGridCell", shapefilename, 1:nshps, 1:ncells,
         xc, yc, dx, dy)
      temp <- tapply(temp$recordArea, temp$cellID, sum)
      extent[as.numeric(names(temp))] <- temp
      prop <- extent/sum(extent)

# Print an error message to indicate unknown shapefile type
   } else {
      stop(paste("\nShapefile type", shp.type, "is not recognized."))
   }

# If a survey design frame object was provided, then  remove the temporary
# shapefile

   if(!is.null(spframe)) {
      file.remove(paste(shapefilename, ".dbf", sep=""), paste(shapefilename,
         ".shp", sep=""), paste(shapefilename, ".shx", sep=""))
   }

# Return results
   list(extent=extent, prop=prop, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
        dx=dx, dy=dy, xc=xc, yc=yc)
}
