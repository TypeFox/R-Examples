simNormMixt <-
function(n, mu.surrog, sd.surrog, pi)
{
  pi<-pi/sum(pi)
  J<-length(mu.surrog)
  r <- rmultinom(1,n,pi)
  unlist(sapply(seq(along=r),function(j) rnorm(r[j], mean=mu.surrog[j], sd = sd.surrog[j])))
}




