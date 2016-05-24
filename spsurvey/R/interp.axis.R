interp.axis <- function(yl.lab, cdfest.l, cdfest.r) {

################################################################################
# Function: interp.axis
# Purpose: Create right side y-axis labels for a CDF plot
# Programmers: Tony Olsen
#              Tom Kincaid
# Date: March 26, 2007
# Last Revised: April 26, 2007
# Description:
#   This function creates right side y-axis labels for a CDF plot.  It assumes
#   that arguments cdfest.l and cdfest.r are strictly increasing.  If argument
#   yl.lab is less than the first cdfest.l vlaue, then the funciton assumes 0
#   for both cdfest.l and cdfest.r.
# Arguments:
#   yl.lab = vector of left side y-axis labels, which are the basis for
#     interpolating cdfest.r values.
#   cdfest.l = vector of CDF estimates corresponding to the left side y-axis.
#   cdfest.r = vector of CDF estimates corresponding to the right side y-axis.
# Output:
#   A numeric vector consisting of the right side y-axis labels.
# Other Functions Required: None
################################################################################

nvec <- 1:length(cdfest.l)
rslt <- numeric(0)
for (j in 1:length(yl.lab)) {
   high <- ifelse(length(nvec[cdfest.l >= yl.lab[j]]) > 0, 
                  min(nvec[cdfest.l >= yl.lab[j]]), nvec[length(nvec)])
   low <- ifelse(length(nvec[cdfest.l <= yl.lab[j]]) >  0,
                 max(nvec[cdfest.l <= yl.lab[j]]), NA)
   if(is.na(low)) {
      ulow <- 0
      plow <- 0
      pdis <- (yl.lab[j] - plow)/(cdfest.l[high] - plow)
    	 rslt[j] <- ulow + pdis * (cdfest.r[high] - ulow)
   } else {
      if( high > low) {
    	    pdis <- (yl.lab[j] - cdfest.l[low])/(cdfest.l[high] - cdfest.l[low])
    	    rslt[j] <- cdfest.r[low] + pdis * (cdfest.r[high] - cdfest.r[low])
      } else {
         rslt[j] <- cdfest.r[high]
      }
   }
}

return(rslt)
}
