rinvgamma <-
function(n, shape, scale) {
  test=rgamma(n, shape=shape, scale=1/scale)
  
  return(1 / rgamma(n, shape=shape, scale=1/scale))
  }
