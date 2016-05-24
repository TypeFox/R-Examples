rdirichlet=function(n, shape){
#Dirichlet generator
  d = length(shape)
  x = rgamma(n*d, rep(shape, n))
  x = matrix(x, nrow = d)
  t(sweep(x, 2, apply(x, 2, sum), "/"))
  }
