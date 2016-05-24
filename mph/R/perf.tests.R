perf.sphere.dim <- function(n=10000, dims = 2:5 ){
  
  dgm = list()
  times = list()
  index = 1
  for(d in dims){
    S = matrix(rnorm( (d+1)*n), nrow=n)
    S = t( scale( t(S), center=FALSE, scale=TRUE) )
# X = noise*matrix(rnorm(D*n), nrow=n)
#   X[,1:(d+1)] = X[,1:(d+1)]  + S
    times[[index ]] = system.time( dgm[[index]] <-  multiscale.rips.ipca(S, -1,
          3 ) ) 
    index = index + 1
  }

  res = list( dgm = dgm, times = times )

}
