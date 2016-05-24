xtable.robustness <-
  function(object, ...){
    
    robustness.object <- object
    

    
    if (class(robustness.object) != "robustness"){
      stop("Wrong object class")
    } 
    
    class(robustness.object) <- "matrix"
    
    if (requireNamespace("xtable", quietly = TRUE)) {
      (print(xtable::xtable(robustness.object)))
    } else {
      stop("xtable library is required for this to function")
    }
    

      
  
    class(robustness.object) <- "robustness"
  }    