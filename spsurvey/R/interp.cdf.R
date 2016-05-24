interp.cdf <- function(pctval, cdfest.p, cdf.value) {

################################################################################
# Function: interp.cdf
# Purpose: Interpolate CDF values at a set of percentiles
# Programmers: Tony Olsen
#              Tom Kincaid
# Date: March 26, 2007
# Last Revised: April 26, 2007
# Description:
#   This function interpolates CDF values at a set of percentiles.  The CDF
#   values can be CDF estimates, CDF confidence bound estimates, or values at
#   which the CDF is estimated (i.e., x-axis values).  It is assumed that
#   arguments cdfest.p and cdf.value are strictly increasing.
# Arguments:
#   pctval = vector of percentiles (expressed as percents) at which the CDF
#     values are to be interpolated.
#   cdfest.p = vector of CDF estimates in terms of proportions.
#   cdf.value = vector of CDF values to be interpolated.
# Output:
#   A numeric vector consisting of the interpolated CDF values.
# Other Functions Required: None
################################################################################

nvec <- 1:length(cdfest.p)
rslt <- numeric(0)
for (j in 1:length(pctval)) {
   high <- ifelse(length(nvec[cdfest.p >= pctval[j]]) > 0, 
                  min(nvec[cdfest.p >= pctval[j]]), NA)
   low <- ifelse(length(nvec[cdfest.p <= pctval[j]]) >  0,
                  max(nvec[cdfest.p <= pctval[j]]), NA)
   if(is.na(high)) {
      rslt[j] <- NA
   } else if(is.na(low)) {
     rslt[j] <- cdf.value[high]
   } else {
      if( high > low) {
    	    pdis <- (pctval[j] - cdfest.p[low])/(cdfest.p[high] - cdfest.p[low])
    	    rslt[j] <- cdf.value[low] + pdis * (cdf.value[high] - cdf.value[low])
      } else {
         rslt[j] <- cdf.value[high]
      }
   }
}

return(rslt)
}

