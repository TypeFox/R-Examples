selectrecordID <- function(rdx, cellID, recordMeasure, recordID, mdm, id) {

################################################################################
# Function: selectrecordID
# Purpose: Select a shapefile record from which to select a sample point
# Programmers: Tom Kincaid
# Date: November 15, 2010
# Description:      
#   This function selects a shapefile record from which a sample point will be
#   selected.  
# Arguments:
#   rdx = the randomized hierarchical address identifying a grid cell that will
#     get a sample point.
#   cellID = the vector of grid cell IDs.
#   recordMeasure = the vector of grid cell shapefile record length for
#     polyline-type shapefiles or shapefile record area for polygon-type
#     shapefiles.
#   recordID = the vector of grid cell shapefile record IDs.
#   mdm = the vector of multidensity multipliers for the shapefile records.
#   id = the vector of shapefile record IDs.
# Results: 
#   The ID of a shapefile record.
# Other Functions Required: None
################################################################################

   ind <- rdx == cellID
   nrec <- sum(ind)
   if(nrec > 1) {
      temp <- recordMeasure[ind] * 
           mdm[match(recordID[ind], id)]
      probs <- temp/sum(temp)
      rslt <- sample(recordID[ind], 1, prob=probs)
   } else {
      rslt <- recordID[ind]
   }
   return(rslt)
}
