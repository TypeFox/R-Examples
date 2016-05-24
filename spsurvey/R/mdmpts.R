mdmpts <- function(mdcaty,n.desired) {

################################################################################
# File: mdmpts.r
# Purpose: Calculate multi-density category multipliers for GRTS applied to
#          points
# Programmer: Tony Olsen
# Date: October 15, 2002
# Input:
#   mdcaty = vector of multi-density category groups for each element in sample
#     frame.
#   n.desired = expected sample size for each category.  Row names must match
#     category names in mdcaty.
# Output:
#   A numeric vector of multipliers that is same length as mdcaty.
################################################################################

   gsum <- table(mdcaty)
   catmatch <- match(names(n.desired),names(gsum),nomatch=0)
   piden <- n.desired/gsum[catmatch]
   mdmpts <- rep(NA,length(mdcaty))
   for(i in names(n.desired))
       mdmpts[mdcaty == i] <- piden[i]

   mdmpts
}


