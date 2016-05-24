DnR <- function(x, d, w, y=rep(-Inf, length(x)) )
{
# inputs should be from  Wdataclean2(), i.e. ordered and weighted.
# y is the truncation times, y do not have to be the same length
# as x, but should be length(y) = sum(w).

allrisk <- cumsumsurv(w)  ## rev(cumsum(rev(w)))  3/2015 MZ
posi <- d == 1
uncenx <- x[posi]
uncenw <- w[posi]
uncenR <- allrisk[posi]

if(any(y > -Inf)) { 
  inde <- function(u, v) { as.numeric(u >= v) }
  uuij <- outer(y, uncenx, FUN="inde")
  trunca <- colSums(uuij)   
  ## as.vector( rowsum( uuij, group= rep(1, length(y))) )  3/2015 MZ
  uncenR <- uncenR - trunca
}

list( times = uncenx, n.risk = uncenR, n.event = uncenw )
}

