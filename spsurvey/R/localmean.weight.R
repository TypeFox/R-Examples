localmean.weight <- function(x, y, prb, nbh=4, vincr=0.00001*abs(mean(y))) {

################################################################################
# Function: localmean.weight
# Programmers: Don Stevens and Tom Kincaid
# Date: September 5, 2001
# Last Revised: September 27, 2013
# Description:
#   This function calculates the index values of neighboring points and
#   associated weights required by the local mean variance estimator.
#   Input:
#      x = x-coordinates for location of the sample points.
#      y = y-coordinates for location of the sample points.
#      prb = inclusion probabilities for the sample points.
#      nbh = number of neighboring points to use in the calculations.
#      vincr = the variance increment for correcting an La.svd error.  The
#        default is 0.00001*abs(mean(y)).
#   Output:
#      An object in list format containing two elements: a matrix named ij 
#      composed of the index values of neighboring points and a vector named gwt
#      composed of weights.
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
   ind <- class(a22) == "try-error"
   iter <- 1
   v <- 0
   while(ind & iter < 11) {
      v <- v + vincr
      xt <- x + rnorm(n, 0, v)
      yt <- y + rnorm(n, 0, v)
      a22 <- localmean.weight2(xt, yt, prb, nbh)
      ind <- class(a22) == "try-error"
      iter <- iter+1
   }
   if(ind)
      stop("\nThe La.svd function terminated with an error.  Add random noise to the \nx-coordinates and y-coordinates.")
   a21 <-  - diag(2/gct) %*% hij %*% a22
   lm <- a21 %*% (1 - smwt)
   gm <- a22 %*% (1 - smwt)
   gwt <- (lm[ij[, 1]] + gm[ij[, 2]])/2 + gwt

# Return the results

   list(ij=ij, gwt=gwt)
}
