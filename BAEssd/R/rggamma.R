rggamma <-
function(n,shape1,rate1,shape2){
  ### Error checking
  # shape2 must be an integer
  shape2 <- as.integer(shape2)
  
  ### Generate y
  y <- rgamma(n,shape=shape1,rate=rate1)
  
  ### Generate X
  out <- rgamma(n,shape=shape2,rate=y)
  
  return(out)
}
