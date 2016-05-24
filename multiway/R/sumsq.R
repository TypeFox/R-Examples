sumsq <- 
  function(X){
    # sum of squares of input object
    
    if(is.list(X)){
      xss <- sum(sapply(X,sumsq))
    } else {
      xss <- sum(X^2)
    }
    
    xss
    
  }