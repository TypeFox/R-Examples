
# hdi method for class stats::density, a list with elements 
#  x = the n coordinates of the points where the density is estimated.
#  y = the estimated density values.

hdi.density <- function(object, credMass=0.95, allowSplit=FALSE, ...) {
  checkCredMass(credMass)
  sorted = sort( object$y , decreasing=TRUE )
  heightIdx = min( which( cumsum( sorted) >= sum(object$y) * credMass ) )
  height = sorted[heightIdx]
  indices = which( object$y >= height )
  # HDImass = sum( object$y[indices] ) / sum(object$y)
  gaps <- which(diff(indices) > 1)
  if(length(gaps) > 0 && !allowSplit) {
    # In this case, return shortest 95% CrI
    warning("The HDI is discontinuous but allowSplit = FALSE;
    the result is a valid CrI but not HDI.")
    cumul <- cumsum(object$y) / sum(object$y)
    upp.poss <- low.poss <- which(cumul < 1 - credMass)
    for (i in low.poss)
      upp.poss[i] <- min(which(cumul > cumul[i] + credMass))
    # all(cumul[upp.poss] - cumul[low.poss] > credMass) # check
    width <- upp.poss - low.poss
    best <- which(width == min(width))  # usually > 1 value due to ties
    result <- c(lower = mean(object$x[low.poss[best]]),
                upper = mean(object$x[upp.poss[best]]))
  } else {
    begs <- indices[c(1, gaps + 1)]
    ends <- indices[c(gaps, length(indices))]
    result <- cbind(begin = object$x[begs], end = object$x[ends])
    if(!allowSplit)  {
      result <- as.vector(result)
      names(result) <- c("lower", "upper")
    }
  }
  attr(result, "credMass") <- credMass
  attr(result, "height") <- height
  return(result)
}
