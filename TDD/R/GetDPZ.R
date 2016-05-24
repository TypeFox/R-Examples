GetDPZ = function(w, dt){
  out = list()
  data('DPZLIST', package = 'TDD', envir = environment())
  DPZLIST = get('DPZLIST', envir = environment())
  if(length(dt) == 1){
    dt = rep(dt, length(w))
  }

  for(i in w){
    if(i > length(DPZLIST)){
      warning(paste(i, 'is not a valid instrument number'))
      next
    }
    L = DPZLIST[[i]]
    wdt = which(sapply(L, function(x)x$dt) == dt[which(w == i)])
    if(length(wdt) == 0){
      warning(paste('No precalculated response for dt =', dt[which(w == i)]))
      next
    }
    out[[which(w == i)]] = L[[wdt]]
  }
  names(out) = names(DPZLIST)[w]
  return(out)
}
