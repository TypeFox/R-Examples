genotypingSavings <-
function(pi.samples, pi.markers) {
  result.joint <- pi.samples + (1-pi.samples)*pi.markers   
  result.savings <- 1 - result.joint
  result.savings
  }
