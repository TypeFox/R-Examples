irsarea <- function (shapefilename=NULL, areaframe, samplesize=100, SiteBegin=1,
   maxtry=1000) {

################################################################################
# Function: irsarea
# Purpose: Select an independent random sample (IRS) of an area resource
# Programmer: Tom Kincaid
# Date: November 30, 2005
# Last Revised: February 23, 2007
# Description:      
#   This function selects an IRS of an area resource.  
# Arguments:
#   shapefilename = name of the input shapefile.  If shapefilename equals NULL,
#     then the shapefile or shapefiles in the working directory are used.  The
#     default is NULL.
#   areaframe = a data frame containing id, mdcaty, area, and mdm.
#   samplesize = number of points to select in the sample.  The default is 100.
#   SiteBegin = first number to start siteID numbering.  The default is 1.
#   maxtry = maximum number of iterations for randomly generating a point within
#     the frame to select a site when type.frame equals "area".  The default
#     is 1000.
# Results: 
#   A data frame of sample points containing: siteID, id, x, y, mdcaty,
#     and weight.
# Other Functions Required:
#   getRecordIDs - C function to obtain the shapefile record IDs for records
#     from which sample points will be selected
#   getShapeBox - C function to obtain the shapefile minimum and maximum values
#     for the x and y coordinates
#   pointInPolygonFile - C function to determine the polygon IDs and the
#     probability values associated with a set of point, where polygons are
#     specified by a shapefile
################################################################################

# Determine IDs for records that will contain sample points

   area.cumsum <- cumsum(areaframe$area*areaframe$mdm)
   samp.pos <- runif(samplesize, 0, area.cumsum[nrow(areaframe)])
   samp.id <- .Call("getRecordIDs", area.cumsum, samp.pos, areaframe$id)

# Pick sample points

   x <- y <- id <- mdm <-numeric(samplesize)
   j <- 1
   for(i in unique(samp.id)) {
      npt <- sum(samp.id == i)
      bp <- !logical(npt)
      ntry <- 0    
      idx <- seq(j, j+npt-1)
      temp <- .Call("getShapeBox", shapefilename, i)
      xmin <- temp$xmin
      ymin <- temp$ymin
      xmax <- temp$xmax
      ymax <- temp$ymax
      while(any(bp) & ntry < maxtry) {
         x[idx[bp]] <- runif(sum(bp), xmin, xmax)
         y[idx[bp]] <- runif(sum(bp), ymin, ymax)
         temp <- .Call("pointInPolygonFile", shapefilename, x[idx[bp]],
            y[idx[bp]], i, areaframe$mdm[areaframe$id == i])
         id[idx[bp]] <- temp$id
         mdm[idx[bp]] <- temp$mdm
         ntry <- ntry + 1
         bp <- mdm[idx] == 0
      }
      j <- j + npt
   }

# When the achieved sample size is less than the desired sample size, remove
# values from the output vectors

   bp <- mdm == 0
   if(sum(!bp) < samplesize) {
      x <- x[!bp]
      y <- y[!bp]
      id <- id[!bp]
      mdm <- mdm[!bp]
   }
   mdcaty <- areaframe$mdcaty[match(id, areaframe$id)]

# Assign Site ID

   siteID <- SiteBegin - 1 + 1:length(x)


# Place Site ID as first column and add weights

   rho <- data.frame(siteID=siteID, id=id, xcoord=x, ycoord=y, mdcaty=mdcaty,
      wgt=1/mdm)
   row.names(rho) <- 1:nrow(rho)

# Return the sample

   rho
}
