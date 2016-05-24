sim.n.chi2 <-
function (nu, lambda) {
  # This function generates a variate from a non-central chi-square distribution 
  # with degree of freedom nu >0 and non-centrality parameter lambda >0 
  
  if( nu > 1){
    print(rchisq(1 , nu-1))
    x <- rchisq(1 , nu-1) + (rnorm(1) + sqrt( lambda ) )^2
  }
  else{
    x <- rchisq(1 , nu + 2*rpois(1, 0.5 * lambda ) )
  }
  
  return( x )
}
