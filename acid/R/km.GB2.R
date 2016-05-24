km.GB2 <-
function(b,a,p,q,k){ #k-th moment function - from Kleiber p.188
  km<- b^k*gamma(p+k/a)*gamma(q-k/a)/(gamma(p)*gamma(q))
  return(km)
}
