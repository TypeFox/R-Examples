soft <-
function(z,lambda){
  # the threshold function of lasso, that is ,the soft threshold operator
  
  ind <- (z>=-lambda) + (z>lambda)
  s <- (z-lambda)*(ind==2) + (z+lambda)*(ind==0)
  return(s)
  
}
