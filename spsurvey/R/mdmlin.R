mdmlin <- function(len, mdcaty, n.desired) {

################################################################################
# File: mdmlin.r
# Purpose: Calculate multi-density category multipliers for GRTS applied to
#          linear networks
# Programmer: Tony Olsen
# Date: April 3, 2003
# Input:
#   len = vector of segment lengths for each segment in sample frame.
#   mdcaty = vector of multi-density category groups for each segment in sample
#     frame.
#   n.desired = expected sample size for each category.  Row names must match
#     category names in mdcaty.
# Output:
#   A numeric vector of multipliers that is same length as len and mdcaty.
################################################################################

   catsum <- tapply(len,mdcaty,sum,na.rm=TRUE)
   catmatch <- match(names(n.desired),names(catsum),nomatch=0)
   piden <- n.desired/catsum[catmatch]
   mdmlin <- rep(NA,length(mdcaty))
   for(i in names(n.desired))
      mdmlin[mdcaty == i] <- piden[i]

   mdmlin
}


