`alowcurve` <- function(p,q,chiquant,n,cc=0.5) {
  y <- (p*q-cc*(1-p*q)/n)^2 - (cc^2*p*q*(p*q-0.5)/n^2 + p^2*q^2*chiquant/n)
  return(y)
}

