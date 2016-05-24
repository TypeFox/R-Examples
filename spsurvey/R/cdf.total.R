cdf.total <- function(z, wgt, val, cluster.ind, cluster, wgt1, popsize) {

################################################################################
# Function: cdf.total
# Programmer: Tom Kincaid
# Date: July 12, 2000
# Last Revised: January 27, 2004
# Description:
#   This function calculates an estimate of the cumulative distribution function
#   (CDF) for the total of a finite or an extensive resource.  The set of
#   values at which the CDF is estimated is supplied to the function.  If the
#   known extent of the resource is provided, the classic ratio estimator is
#   used to calculate the estimate. That estimator is the product of the known
#   extent of the resource and the Horvitz-Thompson ratio estimator, where the
#   latter is the ratio of two Horvitz-Thompson estimators.  The numerator of
#   the ratio estimates the total of the resource equal to or less than a
#   specified value.  The denominator of the ratio estimates the extent of the
#   resource.  If the known extent of the resource is not provided, the Horvitz-
#   Thompson estimator of the total of the resource equal to or less than a
#   specified value is used to calculate the estimate. For a finite resource,
#   size is the number of units in the resource.  For an extensive resource,
#   size is the measure of the resource, i.e., length, area, or volume.  The 
#   function can accomodate single-stage and two-stage samples.
# Arguments:
#   z = the response value for each site.
#   wgt = the final adjusted weight (inverse of the sample inclusion
#     probability) for each site, which is either the weight for a single-
#     stage sample or the stage two weight for a two-stage sample.
#   val = the set of values at which the CDF is estimated.
#   cluster.ind = a logical value that indicates whether the sample is a two-
#     stage sample, where TRUE = a two-stage sample and FALSE = not a two-stage
#     sample.
#   cluster = the stage one sampling unit (primary sampling unit or cluster) 
#     code for each site.
#   wgt1 = the final adjusted stage one weight for each site.
#   popsize = known size of the resource, which is used to perform ratio
#     adjustment to estimators expressed using measurement units for the
#     resource.  For a finite resource, this argument is either the total number
#     of sampling units or the known sum of size-weights.  For an extensive
#     resource, this argument is the measure of the resource, i.e., either known
#     total length for a linear resource or known total area for an areal
#     resource.  For a stratified sample this variable must be a vector
#     containing a value for each stratum and must have the names attribute set
#     to identify the stratum codes.
# Output:
#   The CDF estimate.
# Other Functions Required:
#   None
################################################################################

# Calculate additional required values

   m <- length(val)
   if (cluster.ind) {
      cluster <- factor(cluster)
      ncluster <- length(levels(cluster))
      z.lst <- split(z, cluster)
      wgt2.lst <- split(wgt, cluster)
      wgt1.u <- as.vector(tapply(wgt1, cluster, unique))
   }

# Calculate the cdf estimate

   cdf <- numeric(m)
   if (cluster.ind) {
      for(i in 1:m) {
         temp <- numeric(ncluster)
         for (j in 1:ncluster) {
            temp[j] <- sum(ifelse(z.lst[[j]] <= val[i], wgt2.lst[[j]], 0))
         }
         cdf[i] <- sum(wgt1.u*temp)
      }
   } else {
      for(i in 1:m) {
         cdf[i] <- sum(ifelse(z <= val[i], wgt, 0))
      }
   }

# Adjust the estimate when the size of the resource is known

   if (!is.null(popsize)) {
      if (cluster.ind)
         cdf <- popsize*(cdf/sum(wgt1*wgt))
      else
         cdf <- popsize*(cdf/sum(wgt))
   }

# Return the estimate

   cdf
}
