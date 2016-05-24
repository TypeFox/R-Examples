localmean.var <- function(z, weight.lst) {

################################################################################
# Function: localmean.var
# Programmers: Don Stevens and Tom Kincaid
# Date: October 17, 2000
# Last Revised: September 5, 2001
# Description:
#   This function calculates the local mean variance estimator.
#   Input:
#      z = weighted response values or weighted residual values for the sample
#         points.
#      weight.lst = list from the local mean weight function containing two 
#         elements: a matrix named ij composed of the index values of 
#         neighboring points and a vector named gwt composed of weights.
#   Output: 
#      The local mean estimator of the variance.
#   Other Functions Required: None
################################################################################

# Calculate local means

   zb <- sapply(split(z[weight.lst$ij[, 2]] * weight.lst$gwt, weight.lst$ij[, 1]), sum)

# Calculate the variance estimate

   lmvar <- sum(weight.lst$gwt * (z[weight.lst$ij[, 2]] - zb[weight.lst$ij[, 1]])^2)

# Return the variance estimate

   lmvar
}
