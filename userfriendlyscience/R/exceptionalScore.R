exceptionalScore <- function(x, prob=.025, both=TRUE, silent=FALSE,
                             quantileCorrection = .0001, quantileType = 8) {
  
  belowLower <- x < (quantile(x, probs=min(c(prob, 1-prob)), na.rm=TRUE, type=quantileType) - quantileCorrection);
  aboveUpper <- x > (quantile(x, probs=max(c(prob, 1-prob)), na.rm=TRUE, type=quantileType) + quantileCorrection);
  if (both) {
    return(belowLower | aboveUpper);
  } else {
    if (prob < .5) {
      return(belowLower);
    } else {
      return(aboveUpper);
    }
  }
}
