localmean.df <- function(weight.lst) {

################################################################################
# Function: localmean.df
# Programmer: Tom Kincaid
# Date: April 8, 2003
# Description:
#   This function calculates the degrees of freedom of the local mean variance-
#   covariance estimator.
#   Input:
#      weight.lst = list from the local mean weight function containing two 
#         elements: a matrix named ij composed of the index values of 
#         neighboring points and a vector named gwt composed of weights.
#   Output:
#      The degrees of freedom of the local mean variance-covariance estimator.
#   Other Functions Required: None
################################################################################

# Create the matrix of coefficients used in calculating the local mean variance
# estimator

   n <- max(weight.lst$ij[,1])
   df.mat <- array(0, c(n,n))
   df.mat[weight.lst$ij] <- -weight.lst$gwt
   diag(df.mat) <- 1 + diag(df.mat)

# Calculate the degrees of freedom

   df <- qr(df.mat)$rank

# Return the result

   df
}
