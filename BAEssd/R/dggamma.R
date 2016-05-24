dggamma <-
function(x,shape1,rate1,shape2){
  ### Error checking
  # shape2 must be an integer
  shape2 <- as.integer(shape2)
  
  # positivity
  if(sum(c(x,shape1,rate1,shape2)<=0)>0){
    stop("All parameters and values of x must be positive.")
  }
  
  ### Density
  out <- ((rate1^shape1)/gamma(shape1))*(gamma(shape1+shape2)/gamma(shape2))*
      ((x^(shape2-1))/((rate1+x)^(shape1+shape2)))
  
  return(out)
}
