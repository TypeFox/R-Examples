selectframe <- function(rord, xc, yc, dx, dy, pts) {

################################################################################
# Function: selectframe
# Purpose: Selects all of the points in the frame using the order of cells
#    specified by rord
# Programmer: Tom Kincaid
# Date: February 28, 2006
# Input:
#   rord = the index value for all cells.
#   xc = x-coordinates that define the cells.
#   yc = y-coordinates that define the cells.
#   dx = width of the cells along the x-axis.
#   dy = width of the cells along the y-axis.
#   pts = a data frame containing id values, x-coordinates, y-coordinates, and
#      mdm values.
# Output:
#   The id value for all points in the frame.
################################################################################

   id <- NULL
   for(cel in rord) {
      xr <- c( xc[cel] - dx, xc[cel])
      yr <- c( yc[cel] - dy, yc[cel])
      tstcell <- (xr[1] < pts$x) & (pts$x <= xr[2]) & (yr[1] < pts$y) & (pts$y <= yr[2])
      npt.cell <- length(pts$id[tstcell])
      if(npt.cell == 1) {
         id <- c(id, pts$id[tstcell])
      } else {
         id <- c(id, sample(pts$id[tstcell], npt.cell))
      }
   }
   id
}
