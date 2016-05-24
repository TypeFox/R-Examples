`mean.trapTransect` <-
function(x, ...) {
  return(weighted.mean(x$distances, w=x$seed.count, ...))
}

