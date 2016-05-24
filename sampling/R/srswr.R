"srswr" <-
function(n,N) as.vector(rmultinom(1,n,rep(n/N,times=N)))

