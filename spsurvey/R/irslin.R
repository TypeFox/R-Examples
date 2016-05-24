irslin <- function (shapefilename=NULL, linframe, samplesize=100, SiteBegin=1) {

################################################################################
# Function: irslin
# Purpose: Select an independent random sample (IRS) of a linear resource
# Programmer: Tom Kincaid
# Date: November 17, 2005
# Last Revised: March 7, 2007
# Description:      
#   This function selects an IRS of a linear resource.  
# Arguments:
#   shapefilename = name of the input shapefile.  If shapefilename equals NULL,
#     then the shapefile or shapefiles in the working directory are used.  The
#     default is NULL.
#   linframe = a data frame containing id, mdcaty, len, and mdm.
#   samplesize = number of points to select in the sample.  The default is 100.
#   SiteBegin = first number to start siteID numbering.  The default is 1.
# Results: 
#   A data frame of sample points containing: siteID, id, x, y, mdcaty,
#   and weight.
# Other Functions Required: None
################################################################################

# Pick sample points

   len.cumsum <- cumsum(linframe$len*linframe$mdm)
   samp.pos <- runif(samplesize, 0, len.cumsum[nrow(linframe)])
   ordr <- rank(samp.pos)
   samp.pos <- sort(samp.pos)
   temp <- .Call("linSampleIRS", shapefilename, len.cumsum, samp.pos,
      linframe$id, linframe$len, linframe$mdm)
   temp$id <- temp$id[ordr]
   temp$x <- temp$x[ordr]
   temp$y <- temp$y[ordr]
   mdcaty <- linframe$mdcaty[match(temp$id, linframe$id)]
   mdm <- linframe$mdm[match(temp$id, linframe$id)]

# Assign Site ID

   siteID <- SiteBegin - 1 + 1:length(temp$id)

# Create the output data frame

   rho <- data.frame(siteID=siteID, id=temp$id, xcoord=temp$x, ycoord=temp$y,
      mdcaty=mdcaty, wgt=1/mdm)
   row.names(rho) <- 1:nrow(rho)

# Return the sample

   rho
}
