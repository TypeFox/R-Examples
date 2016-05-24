mdmarea <- function(area, mdcaty, n.desired) {

################################################################################
# File: mdmarea.r
# Purpose: Calculate multi-density category multipliers for GRTS applied to
#          polygonal areas
# Programmer: Tom Kincaid
# Date: October 12, 2004
# Input:
#   area = vector of polygon areas for each polygon in the sample frame.
#   mdcaty = vector of multi-density category names for each polygon in the
#     sample frame.
#   n.desired = expected sample size for each category.  Row names must match
#     category names in mdcaty.
# Output:
#   A numeric vector of multipliers that is same length as area and mdcaty.
################################################################################

   catsum <- tapply(area, mdcaty, sum, na.rm=TRUE)
   catmatch <- match(names(n.desired),names(catsum),nomatch=0)
   piden <- n.desired/catsum[catmatch]
   mdmarea <- rep(NA, length(mdcaty))
   for(i in names(n.desired))
      mdmarea[mdcaty == i] <- piden[i]

   mdmarea
}


