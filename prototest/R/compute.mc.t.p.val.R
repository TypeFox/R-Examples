compute.mc.t.p.val <-
function(ts, limits){
  # limits
  top.large = min (ts, limits[2])
  top.small = max(-ts, limits[1])
  
  # probabilities
  top.p.large = pnorm (top.large, 0, 1, log.p=TRUE)
  top.p.small = pnorm (top.small, 0, 1, log.p=TRUE)
  bottom.p.large = pnorm (limits[2], 0, 1, log.p=TRUE)
  bottom.p.small = pnorm (limits[1], 0, 1, log.p=TRUE)
  
  # probability
  numer = top.p.large + log (1 - exp(top.p.small-top.p.large))
  denom = bottom.p.large + log (1 - exp(bottom.p.small-bottom.p.large))
  
  return (1-exp(numer-denom))
}
