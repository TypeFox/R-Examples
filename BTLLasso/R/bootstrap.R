boot.BTLLasso <- function(model, B=500, lambda = NULL, cores = 1, trace = TRUE, 
                          trace.cv = TRUE, control = BTLLasso.ctrl()){
  
  Y <- model$Y
  X <- model$X
  I <- ncol(Y)
  m <- (1 + sqrt(1+8*I))/2
  p <- ncol(X)
  n <- nrow(Y)
  n.theta <- model$n.theta
  
  folds <- model$folds
  
  if(is.null(lambda)){
    lambda <- model$lambda
  }
  
  
  boot.fun <- function(b){
      cat("Bootstrap sample:",b,"out of", B,"\n")
#       options(error = recover) 
#       options(showErrorCalls = T) 

      fullrank <- FALSE
      
      while(!fullrank){
      sample.b <- sample(1:n,replace=TRUE)
      
      X.b <- X[sample.b,]
      Y.b <- Y[sample.b,]
      
      fullrank <- rankMatrix(X.b)[[1]] == ncol(X.b)
      }
      
      model.b <- try(cv.BTLLasso(Y=Y.b,X=X.b,folds=folds,lambda=lambda,cores=1,trace=trace,
                                trace.cv = trace.cv, control = control))
      if(inherits(model.b, "try-error")){
        coef.b <- NA
        lambda.b <- NA
      }else{
      coef.b <- model.b$coefs[which.min(model.b$deviances),]
      lambda.b <- lambda[which.min(model.b$deviances)]
      }
    return(list(coef.b = coef.b, lambda.b = lambda.b))
  }
  

  if(cores>1){
    cl <- makeCluster(cores,outfile="")

    clusterExport(cl, varlist = c("X","Y","folds","lambda", "B","n"), 
                  envir = sys.frame(sys.nframe()))
    outputB <- parLapply(cl, seq(B), boot.fun)
    stopCluster(cl)
  }else{
    outputB <- lapply(seq(B), boot.fun)
  }

  
  
  estimatesB <- matrix(0,ncol=(m-1)*(p+1)+n.theta,nrow = B)
  lambdaB <- c()
  for(b in 1:B){    
    if(any(is.na(outputB[[b]]$coef.b))){
      estimatesB[b,] <- rep(NA,ncol(estimatesB))
      lambdaB[b] <- NA
    }else{
      estimatesB[b,] <- outputB[[b]]$coef.b
      lambdaB[b] <- outputB[[b]]$lambda.b
    }
  }
  
  estimatesBrepar <- expand.coefs(estimatesB,m=model$m, n.theta = model$n.theta)
  
  
  conf.ints <- apply(estimatesB,2,quantile, probs = c(0.025,0.975), type=1)
  conf.ints.repar <- apply(estimatesBrepar,2,quantile, probs = c(0.025,0.975), type=1)
  
  lambda.min.alert <- any(lambdaB==min(lambda))
  lambda.max.alert <- any(lambdaB==max(lambda))
  
returns <- list(cv.model = model, estimatesB = estimatesB, estimatesBrepar = estimatesBrepar, lambdaB = lambdaB, 
                conf.ints = conf.ints, conf.ints.repar = conf.ints.repar, 
                lambda.max.alert = lambda.max.alert, lambda.min.alert = lambda.min.alert)

class(returns) <- "boot.BTLLasso"

  return(returns)
  
}