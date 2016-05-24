print.summary.bigssa <- 
  function(x,digits=4,...){
    
    x$info <- round(x$info,digits)
    if(is.null(x$residuals)){
      cat("\nCall:\n")
      print(x$call)
      cat("\nError Std Dev Estimate:\n")
      cat(x$sigma," on ",as.numeric(x$n-x$df)," degrees of freedom","\n")
      cat("\nFit Statistics:")
      cat("\nGCV:  ",x$info[1])
      cat("\nR^2:  ",x$info[2])
      cat("\nAIC:  ",x$info[3])
      cat("\nBIC:  ",x$info[4],"\n ")
      cat("\nSmoothing Parameters:\n")
      print(c(lambda=x$lambda,x$gammas))
      cat("\n")
    } else {
      ehat <- round(quantile(x$residuals),digits)
      names(ehat) <- c("Min","1Q","Median","3Q","Max")
      cat("\nCall:\n")
      print(x$call)
      cat("\nResiduals:\n")
      print(ehat)
      cat("\nError Std Dev Estimate:\n")
      cat(x$sigma," on ",as.numeric(x$n-x$df)," degrees of freedom","\n")
      cat("\nFit Statistics:")
      cat("\nGCV:  ",x$info[1])
      cat("\nR^2:  ",x$info[2])
      cat("\nAIC:  ",x$info[3])
      cat("\nBIC:  ",x$info[4],"\n ")
      cat("\nSmoothing Parameters:\n")
      print(c(lambda=x$lambda,x$gammas))
      cat("\n")
    }
    
  }