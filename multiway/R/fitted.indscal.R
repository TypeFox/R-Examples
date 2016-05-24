fitted.indscal <-
  function(object,...){
    # Calculates Fitted Values (arrays) for fit INDSCAL Models
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: April 9, 2015
    
    mydim <- c(rep(nrow(object$B),2),nrow(object$C))
    array(tcrossprod(object$B,krprod(object$C,object$B)),dim=mydim)
    
  }