.weighted.sum = function(x, w, ..., na.rm=FALSE) {
  if(isTRUE(na.rm)) {
    ind = !is.na(x)
    x = x[ind]
    w = w[ind]
  }
  return(sum(x*w))
}

.globalFitness = function(fitness, opt) {
  # must receive fitness and optionally weight argument
  aggFn = match.fun(opt$aggFn)
  out   = apply(fitness, 1, aggFn, w=opt$weights)
  return(out)

}

