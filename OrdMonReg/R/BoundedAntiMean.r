BoundedAntiMean <- function(y, w, a = NA, b = NA){
#
# Pool-adjacent-violaters-algorithm for a antitonic weightedMA
# mean, bounded by a lower bound a and an upper bound b.
#
# Input:
#   - y : data points
#   - w : corresponding weights
#   - a : vector of lower bound
#   - b : vector of upper bound
#
# Fadoua Balabdaoui and Kaspar Rufibach, October 2008

res <- -BoundedIsoMean(-y, w, a = -b, b = -a)
return(res)

}
