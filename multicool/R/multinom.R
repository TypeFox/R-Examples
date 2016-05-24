multinom = function(x, counts = FALSE){
  
  u = NULL
  
  if(!counts){
    u = as.vector(table(x))
  }else{
    ## make sure x is a vector of counts
    is.wholenumber = function(x, tol = .Machine$double.eps^0.5){abs(x - round(x)) < tol}
    
    if(any(!is.wholenumber(x)) || any(x < 0))
      stop("if counts == TRUE then all elements of x must be integer and >= 0")
    
    u = x
  }
  

  r = .Call('multicool_multinomCoeff', PACKAGE = 'multicool', u)
  
  return(r)
}