localmean.weight2 <- function(x, y, prb, nbh) {

################################################################################
# Function: localmean.weight2
# Programmers: Don Stevens and Tom Kincaid
# Date: September 17, 2013
# Description:
#   This function calculates the initial section of the localmean.weight
#   function and serves to allow recovery from an error in the singular value
#   decomposition function (La.svd) that is called by the generalized inverse
#   function (ginv) in the MASS package.
#   Input:
#      x = x-coordinates for location of the sample points.
#      y = y-coordinates for location of the sample points.
#      prb = inclusion probabilities for the sample points.
#      nbh = number of neighboring points to use in the calculations.
#   Output:
#      Either an object of class "try-error" when the ginv function terminates
#      with an error or a generalized inverse matrix when the ginv function
#      terminates normally.
################################################################################

   n <- length(x)

# Calculate indices of nearest neighbors

   dst <- as.matrix(dist(cbind(x, y), diag = TRUE, upper = TRUE))
   idx <- apply(dst, 2, order)[1:nbh,]

# Make neighbors symmetric

   jdx <- rep(1:n, rep(nbh, n))
   kdx <- unique(c((jdx - 1) * n + idx, (idx - 1) * n + jdx)) - 1
   ij <- cbind((kdx) %/% n + 1, (kdx) %% n + 1)
   ij <- ij[order(ij[, 1], dst[ij]),]

# Apply linear taper to the  inverse probability weights

   gct <- tabulate(ij[, 1])
   gwt <- numeric(0)
   for(i in 1:n)
      gwt <- c(gwt, 1 - (1:gct[i] - 1)/(gct[i]))
   gwt <- gwt/prb[ij[, 2]]

# Normalize to make true average

   smwt <- sapply(split(gwt, ij[, 1]), sum)
   gwt <- gwt/smwt[ij[, 1]]
   smwt <- sapply(split(gwt, ij[, 2]), sum)

# Make weights doubly stochastic

   hij <- matrix(0, n, n)
   hij[ij] <- 0.5
   a22 <- try(ginv(diag(gct/2) - hij %*% diag(2/gct) %*% hij), TRUE)

# Return the result

   a22
}
