irspts <- function(ptsframe, samplesize=100, SiteBegin=1) {

################################################################################
# Function: irspts
# Purpose: Select an independent random sample (IRS) of a finite resource
# Programmer: Tom Kincaid
# Date: November 16, 2005
# Last Revised: December 15, 2005
# Description:
#   This function selects an IRS of a finite resource (discrete points).  
# Arguments:
#   ptsframe = a data frame containing id, x, y, mdcaty, and mdm.
#   samplesize = number of points to select in the sample.  The default is 100.
#   SiteBegin = first number to start siteID numbering.  The default is 1.
# Results: 
#   A data frame of sample points containing: siteID, id, x, y, mdcaty,
#   and weight.
# Other Functions Required: None
################################################################################

# Pick sample points

   if(nrow(ptsframe) <= samplesize) {
      id <- ptsframe$id
   } else {
      id <- sample(ptsframe$id, samplesize, prob=ptsframe$mdm)
   }
   temp <- ptsframe[match(id, ptsframe$id), ]
   temp$id <- factor(temp$id)

# Assign Site ID

   siteID <- SiteBegin - 1 + 1:nrow(temp)

# Create the output data frame

   rho <- data.frame(siteID=siteID, id=temp$id, xcoord=temp$x, ycoord=temp$y,
      mdcaty=temp$mdcaty, wgt=1/temp$mdm)
   row.names(rho) <- 1:nrow(rho)

# Return the sample

   rho
}
