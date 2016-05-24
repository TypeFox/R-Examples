localmean.cov <- function(zmat, weight.lst) {

################################################################################
# Function: localmean.cov
# Programmer: Tom Kincaid
# Date: November 2, 2000
# Last Revised: September 5, 2001
# Description:
#   This function calculates the variance-covariance matrix using the local mean
#   estimator.
#   Input:
#      zmat = matrix of weighted response values or weighted residual values for
#         the sample points.
#      weight.lst = list from the local mean weight function containing two 
#         elements: a matrix named ij composed of the index values of 
#         neighboring points and a vector named gwt composed of weights.
#   Output: 
#      The local mean estimator of the variance-covariance matrix.
#   Other Functions Required: None
################################################################################

# Calculate additional required values

   temp <- dim(zmat)
   m <- temp[2]

# Initialize the results matrix

   lmvar <- array(0, c(m, m))

# Begin loops for variance/covariance calculations

   for (k in 1:m) {

      for (l in k:m) {

         z1 <- zmat[, k]
         z2 <- zmat[, l]

# Calculate local means

         zb1 <- sapply(split(z1[weight.lst$ij[, 2]] * weight.lst$gwt, weight.lst$ij[, 1]), sum)
         zb2 <- sapply(split(z2[weight.lst$ij[, 2]] * weight.lst$gwt, weight.lst$ij[, 1]), sum)

# Calculate the variance or covariance estimate

         lmvar[k, l] <- sum(weight.lst$gwt * (z1[weight.lst$ij[, 2]] - zb1[weight.lst$ij[, 1]]) * (z2[weight.lst$ij[, 2]] - zb2[weight.lst$ij[, 1]]))

      }

# Assign estimates that already have been calculated

      if (k > 1) {
         lmvar[k, 1:(k-1)] <- lmvar[1:(k-1), k]
      }
   }

# Return the variance/covariance estimate

   lmvar
}
