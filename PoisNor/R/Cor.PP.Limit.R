Cor.PP.Limit <-
function(lamvec){
  samples=100000
  u = runif(samples, 0, 1)
  lambda1=lamvec[1]
  lambda2=lamvec[2]
  maxcor=cor(qpois(u, lambda1), qpois(u, lambda2))
  mincor=cor(qpois(u, lambda1), qpois(1-u, lambda2))
  return( c(mincor,maxcor) )
}
