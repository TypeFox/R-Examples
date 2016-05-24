dsample.strata <-
function(sid,Ns,ns)
{ ## sid : vector strata ids
  ## Ns  : mos each strata
  ## ns  : sample sizes each strata
  ## value : id sampled from each strata
  N <- length(sid)
  oid <- sid + runif(N)
  oo <- order(oid)
  i0 <- c(0, cumsum(Ns)[-length(Ns)])
  i0 <- rep(i0,ns)
  i1 <- sequence(ns)
  ii <- i0 + i1
  oo[ii]
}
