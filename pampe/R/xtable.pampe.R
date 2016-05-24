xtable.pampe <-
  function(object, ttype, ...){
    
    pampe.object <- object
    
    if (class(pampe.object) != "pampe"){
      stop("Wrong object class")
    } 
    
    if (ttype != "model" & ttype != "treatment"){
      stop("Wrong table type")
    } 
    
    if (requireNamespace("xtable", quietly = TRUE)) {

    
   
    if (ttype=="model"){
    
      print(xtable::xtable(pampe.object$model))
      
    }
    
    if (ttype=="treatment"){
      apt.table <- cbind(pampe.object$counterfactual,
                         pampe.object$counterfactual[,1] - pampe.object$counterfactual[,2])
      ##Give names
      colnames(apt.table) <- c("Actual", "Predicted", "Treated")
    
      print(xtable::xtable(apt.table))
    }
    
    }
    
    else {
      stop("xtable library is required for this to function")
    }
    
    
  }    