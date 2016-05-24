plotresiduals <-
function(X,Y,betaresiduals,type){
  
  if(is.null(betaresiduals)){
   stop("Residual data is not included")
  }
  
  if(is.null(X)|is.null(Y)){
   stop("There is no data")
  }
  
  if(max(Y)> 1 | min(Y) < 0){
    stop("The data is not in the (0,1) interval")
  }
  
  if(!is.matrix(X)){
   stop("The variables must be a matrix")
  }
  
  
  Y <- as.matrix(Y)
  
  
  if(type==1 | type==5){
    
    plot(betaresiduals$deviance, main="Deviance residuals")
    
    for (i in 1:(ncol(X)-1)){
      plot(betaresiduals$deviance,X[,i+1], main=paste("Deviance residuals vs X", i),xlab=paste("X", i), ylab=paste("Deviance Residuals"))
    }
    for(i in 1:ncol(Y)){
      plot(betaresiduals$deviance, Y[,i], main=paste("Deviance Residuals vs Y",i), xlab=paste("Y", i), ylab=paste("Deviance Residuals"))
    }
    
  } 
  else if(type==2 | type==5){
    
    plot(betaresiduals$pearson, main="Pearson residuals")
    
    for (i in 1:(ncol(X)-1)){
      plot(betaresiduals$pearson,X[,i+1], main=paste("Pearson residuals vs X", i),xlab=paste("X", i), ylab=paste("Pearson Residuals"))
    }
    for(i in 1:ncol(Y)){
      plot(betaresiduals$pearson, Y[,i], main=paste("Pearson Residuals vs Y",i), xlab=paste("Y", i), ylab=paste("Pearson Residuals"))
    }
  }
   else if(type==3| type==5){
    
    plot(betaresiduals$std.pearson, main="Standardized Pearson residuals")
    
    for (i in 1:(ncol(X)-1)){
      plot(betaresiduals$std.pearson,X[,i+1], main=paste("Standardized Pearson residuals vs X", i),xlab=paste("X", i), ylab=paste("Pearson Residuals"))
    }
    for(i in 1:ncol(Y)){
      plot(betaresiduals$std.pearson, Y[,i], main=paste("Standardized Pearson residuals vs Y",i), xlab=paste("Y", i), ylab=paste("Pearson Residuals"))
    }
  }
  else {
    
    plot(betaresiduals$abs, main="Residuals")
    
    for (i in 1:(ncol(X)-1)){
      plot(betaresiduals$abs,X[,i+1], main=paste("Residuals vs X", i),xlab=paste("X", i), ylab=paste("Pearson Residuals"))
    }
    for(i in 1:ncol(Y)){
      plot(betaresiduals$abs, Y[,i], main=paste("Residuals vs Y",i), xlab=paste("Y", i), ylab=paste("Pearson Residuals"))
    }
  }
}
