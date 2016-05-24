sbcsamp <- function(sp.sample, sbc.frame=NULL, dx=NULL, dy=NULL,
   xc=NULL, yc=NULL) {

################################################################################
# Function: sbcsamp
# Purpose: Calculate spatial balance grid cell extent and proportions for a
#          survey design
# Programmer: Tom Kincaid
# Date: September 29, 2011
# Description:      
#   This function calculates spatial balance grid cell extent and proportions
#   for a survey design.  The user must provide either sbc.frame or values for
#   dx, dy, xc, and yc.
# Arguments:
#   sp.sample = the sp package object of class "SpatialPointsDataFrame" produced by
#     the grts or irs functions that contains survey design information.
#   sbc.frame = the object created by the sbcframe function.  The default is
#     NULL.
#   dx = grid cell x-coordinate increment value.  The default is NULL.
#   dy = grid cell y-coordinate increment value.  The default is NULL.
#   xc = vector of grid cell x-coordinates.  The default is NULL.
#   yc = vector of grid cell y-coordinates.  The default is NULL.
# Results: 
#   A list containing the following components: (1) extent - the sample extent
#   for each grid cell and (2) prop - the sample proportion for each grid cell
# Other Functions Required:
#   readShapeFilePts - C function to read the shp file of a point shapefile and
#     return a data frame containing the x-coordinates and y-coordinates for
#     elements in the frame
#   cell.wt - calculates number of points in a cell for a points object
################################################################################

# Obtain the sample x-coordinates and y-coordinates from the sp.sample object
   xcoord <- sp.sample@data$xcoord
   ycoord <- sp.sample@data$ycoord

# If the sbc.frame object was provided, obtain values for dx, dy, xc, and yc
   if(!is.null(sbc.frame)) {
      dx <- sbc.frame$dx
      dy <- sbc.frame$dy
      xc <- sbc.frame$xc
      yc <- sbc.frame$yc
   }

# Calculate grid cell extent and proportion
   ncells <- length(xc)
   ptsframe <- data.frame(x=xcoord, y=ycoord, mdm=1)
   extent <- sapply(1:ncells, cell.wt, xc, yc, dx, dy, ptsframe)
   prop <- (extent/sum(extent))

# Return results
   list(extent=extent, prop=prop)
}
