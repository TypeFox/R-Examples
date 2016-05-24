print.summary.pampe <-
  function(x, ... ){
    
    pampe.object <- x
    
    if (class(pampe.object) != "summary.pampe"){
      stop("Wrong object class")
    } 
    
    time.tr <- pampe.object$time.tr
    
    cat("Selected controls:\n ")
    
    cat(paste(pampe.object$controls[-length(pampe.object$controls)],",", sep=""))
    cat(paste(" and ", pampe.object$controls[length(pampe.object$controls)],".\n",sep=""))
    
    cat("\n \nTime-average estimated treatment effect:\n ")
    
    cat(pampe.object$avg.tr.effect)
    
    cat("\n\nOptimal model estimation results:\n\n")
    
    printCoefmat(pampe.object$model$coefficients)
    
    cat("\nResidual standard error: ")
    cat(pampe.object$model$sigma)
    cat(" on ")
    cat(pampe.object$model$df)
    cat(" degrees of freedom")
    cat("\nMultiple R-squared: ")
    cat(pampe.object$model$r.squared)
    cat(",     Adjusted R-squared: ")
    cat(pampe.object$model$adj.r.squared)
    cat("\nF-statistic: ")
    cat(round(pampe.object$model$fstatistic[1],2))
    cat(" on ")
    cat(round(pampe.object$model$fstatistic[2]))
    cat(" and ")
    cat(round(pampe.object$model$fstatistic[3]))
    cat(" DF, p-value: ")
    cat(pf(pampe.object$model$fstatistic[1], pampe.object$model$fstatistic[2], pampe.object$model$fstatistic[3], lower.tail=F))
    
   
    
  }



