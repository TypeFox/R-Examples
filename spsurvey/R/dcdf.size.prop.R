dcdf.size.prop <- function(g, wgt, cluster.ind, cluster, wgt1, swgt, swgt1) {

################################################################################
# Function: dcdf.size.prop
# Programmer: Tom Kincaid
# Date: December 3, 2002
# Last Revised: June 4, 2008
# Description:
#   This function calculates an estimate of the size-weighted, deconvoluted  
#   cumulative distribution function (CDF) for the proportion of a discrete
#   resource. The simulation extrapolation deconvolution method (Stefanski and
#   Bay, 1996) is use to deconvolute measurement error variance from the 
#   response.  The Horvitz-Thompson ratio estimator, i.e., the ratio of two  
#   Horvitz-Thompson estimators, is used to calculate the estimate.  The 
#   numerator of the ratio estimates the size-weighted total of the resource 
#   equal to or less than a specified value.  The denominator of the ratio 
#   estimates the sum of the size-weights for the resource.  The function can 
#   accomodate single-stage and two-stage samples.
# Arguments:
#   g = values of the deconvolution function g(.) evaluated at a specified value
#     for the response value for each site.
#   wgt = the final adjusted weight (inverse of the sample inclusion
#     probability) for each site, which is either the weight for a single-
#     stage sample or the stage two weight for a two-stage sample.
#   cluster.ind = a logical value that indicates whether the sample is a two-
#     stage sample, where TRUE = a two-stage sample and FALSE = not a two-stage
#     sample.
#   cluster = the stage one sampling unit (primary sampling unit or cluster) 
#     code for each site.
#   wgt1 = the final adjusted stage one weight for each site.
#   swgt = the size-weight for each site, which is the stage two size-weight for
#     a two-stage sample.
#   swgt1 = the stage one size-weight for each site.
# Output:
#   The deconvoluted CDF estimate.
# Other Functions Required:
#   None
################################################################################

# Calculate additional required values

   wgt <- wgt*swgt
   if (cluster.ind) {
      cluster <- factor(cluster)
      ncluster <- length(levels(cluster))
      wgt2.lst <- split(wgt, cluster)
      wgt1 <- wgt1*swgt1
      wgt1.u <- as.vector(tapply(wgt1, cluster, unique))
   }

# Calculate the cdf estimate

   if (cluster.ind) {
      temp <- array(0, c(ncluster, dim(g[[1]])[2]))
      for (i in 1:ncluster) {
         temp[i,] <- apply(g[[i]]*wgt2.lst[[i]], 2, sum)
      }
      cdf <- apply(wgt1.u*temp, 2, sum)
   } else {
      cdf <- apply(wgt*g, 2, sum)
   }

   if (cluster.ind)
      cdf <- cdf/sum(wgt1*wgt)
   else
      cdf <- cdf/sum(wgt)

# Return the estimate

   cdf
}
