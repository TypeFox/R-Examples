qFun = function(probs, x, px){
  Px = cumsum(px)
  
  r = rep(0, length(probs))
  r[probs < 0 | probs > 1] = NA
  r[probs = 0] = min(x)
  r[probs = 1] = max(x)
  
  qp = function(p){
    d = which(Px - p >= 0)[1]
    x[d]
  }
  
  i = probs > 0 & probs < 1
  r[i] = sapply(probs[i], qp)
  return(r)
}