cdf.test.prop <- function(z, wgt, bounds, cluster.ind, cluster, wgt1) {

################################################################################
# Function: cdf.test.prop
# Programmer: Tom Kincaid
# Date: November 5, 2007
# Description:
#   This function calculates a size-weighted estimate of the population 
#   proportions in a set of intervals (classes).  The set of values defining the 
#   upper bound of each class is supplied to the function.  The Horvitz-Thompson 
#   ratio estimator, i.e., the ratio of two Horvitz-Thompson estimators, is used 
#   to calculate the estimate.  The numerator of the ratio estimates the total 
#   of the resource within a class.  The denominator of the ratio estimates the 
#   size of the resource.  For a finite resource size is the number of units in 
#   the resource.  For an extensive resource size is the extent (measure) of the 
#   resource, i.e., length, area, or volume.  The function can accomodate single 
#   stage and two-stage samples.  
# Arguments:
#   z = the response value for each site.
#   wgt = the final adjusted weight (inverse of the sample inclusion
#     probability) for each site, which is either the weight for a single-stage
#     sample or the stage two weight for a two-stage sample.
#   bounds = upper bounds for calculating classes for the CDF.
#   cluster.ind = a logical value that indicates whether the sample is a two-
#     stage sample, where TRUE = a two-stage sample and FALSE = not a two-stage
#     sample.
#   cluster = the stage one sampling unit (primary sampling unit or cluster) 
#     code for each site.
#   wgt1 = the final adjusted stage one weight for each site.
# Results:
#   The class proportions estimate.
# Other Functions Required: None
################################################################################

# Calculate additional required values

   m <- length(bounds)
   ubound <- rep(bounds, m)
   lbound <- c(-1e10, bounds[-m])
   if(cluster.ind) {
      cluster <- factor(cluster)
      ncluster <- length(levels(cluster))
      z.lst <- split(z, cluster)
      wgt2.lst <- split(wgt, cluster)
      wgt1.u <- as.vector(tapply(wgt1, cluster, unique))
   }

# Calculate the class proportions estimate

   phat <- numeric(m)
   if(cluster.ind) {
      for(i in 1:m) {
         temp <- numeric(ncluster)
         for(j in 1:ncluster) {
            temp[j] <- sum(ifelse(z.lst[[j]] <= ubound[i], wgt2.lst[[j]], 0) -
               ifelse(z.lst[[j]] <= lbound[i], wgt2.lst[[j]], 0))
         }
         phat[i] <- sum(wgt1.u*temp)
      }
   } else {
      for(i in 1:m) {
         phat[i] <- sum(ifelse(z <= ubound[i], wgt, 0) - ifelse(z <= lbound[i],
            wgt, 0))
      }
   }

   if(cluster.ind)
      phat <- phat/sum(wgt1*wgt)
   else
      phat <- phat/sum(wgt)

# Return the estimate

   phat
}
