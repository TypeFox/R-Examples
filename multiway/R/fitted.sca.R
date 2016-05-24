fitted.sca <-
  function(object,...){
    # Calculates Fitted Values (lists) for fit SCA Models
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: April 9, 2015
    
    sapply(object$D,function(d,b){tcrossprod(d,b)},b=object$B)
    
  }