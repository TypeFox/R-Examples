cell.wt <- function(cel, xc, yc, dx, dy, pts) {

################################################################################
# Function: cell.wt
# Purpose: Calculates the total inclusion probability for a cell.
# Programmer: Tony Olsen
# Date: October 27, 2004
# Input:
#   cel = the index value for a cell.
#   xc = x-coordinates that define the cells.
#   yc = y-coordinates that define the cells.
#   dx = width of the cells along the x-axis.
#   dy = width of the cells along the y-axis.
#   pts = a data frame containing x-coordinates, y-coordinates, and mdm values.
# Output:
#   The total inclusion probability for a cell.
################################################################################

   xr <- c( xc[cel] - dx, xc[cel])
   yr <- c( yc[cel] - dy, yc[cel])
   tstcell <- (xr[1] < pts$x) & (pts$x <= xr[2]) & (yr[1] < pts$y) & (pts$y <= yr[2])
   wt <- sum(pts$mdm[tstcell])
   wt
}
