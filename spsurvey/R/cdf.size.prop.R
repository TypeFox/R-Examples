cdf.size.prop <- function(z, wgt, val, cluster.ind, cluster, wgt1, swgt,
   swgt1) {

################################################################################
# Function: cdf.size.prop
# Programmer: Tom Kincaid
# Date: July 12, 2000
# Last Revised: June 3, 2008
# Description:
#   This function calculates an estimate of the size-weighted cumulative
#   distribution function (CDF) for the proportion of a finite resource. The 
#   set of values at which the CDF is estimated is supplied to the function. The
#   Horvitz-Thompson ratio estimator, i.e., the ratio of two Horvitz-Thompson 
#   estimators, is used to calculate the estimate.  The numerator of the ratio
#   estimates the size-weighted total of the resource equal to or less than a
#   specified value.  The denominator of the ratio estimates the sum of the
#   size-weights for the resource.  The function can accomodate single-stage and 
#   two-stage samples.
# Arguments:
#   z = the response value for each site.
#   wgt = the final adjusted weight (inverse of the sample inclusion
#     probability) for each site, which is either the weight for a single-
#     stage sample or the stage two weight for a two-stage sample.
#   swgt = the size-weight for each site, which is the stage two size-weight for
#     a two-stage sample.
#   val = the set of values at which the CDF is estimated.
#   cluster.ind = a logical value that indicates whether the sample is a two-
#     stage sample, where TRUE = a two-stage sample and FALSE = not a two-stage
#     sample.
#   cluster = the stage one sampling unit (primary sampling unit or cluster) 
#     code for each site.
#   wgt1 = the final adjusted stage one weight for each site.
#   swgt1 = the stage one size-weight for each site.
# Output:
#   The size-weighted CDF estimate.
# Other Functions Required:
#   None
################################################################################

# Calculate additional required values

   m <- length(val)
   wgt <- wgt*swgt
   if(cluster.ind) {
      cluster <- factor(cluster)
      ncluster <- length(levels(cluster))
      z.lst <- split(z, cluster)
      wgt2.lst <- split(wgt, cluster)
      wgt1 <- wgt1*swgt1
      wgt1.u <- as.vector(tapply(wgt1, cluster, unique))
   }

# Calculate the cdf estimate

   cdf <- numeric(m)
   if(cluster.ind) {
      for(i in 1:m) {
         temp <- numeric(ncluster)
         for(j in 1:ncluster) {
            temp[j] <- sum(ifelse(z.lst[[j]] <= val[i], wgt2.lst[[j]], 0))
         }
         cdf[i] <- sum(wgt1.u*temp)
      }
   } else {
      for(i in 1:m) {
         cdf[i] <- sum(ifelse(z <= val[i], wgt, 0))
      }
   }

   if(cluster.ind)
      cdf <- cdf/sum(wgt1*wgt)
   else
      cdf <- cdf/sum(wgt)

# Return the estimate

   cdf
}
