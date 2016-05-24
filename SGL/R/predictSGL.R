predictSGL = function(x,newX,lam){
  cvobj = x

  X <- newX

  if(!is.null(x$X.transform)){
    if(is.matrix(X)){
      X <- t(t(newX) - x$X.transform$X.means)
      X <- t(t(X) / x$X.transform$X.scale)
    }
    if(is.vector(X)){
      X <- X - x$X.transform$X.means
      X <- X / x$X.transform$X.scale
    }  
  }

  intercept <- 0

  if(!is.null(x$intercept)){
    if(x$type == "linear"){
      intercept <- x$intercept
    }
    if(x$type == "logit"){
      intercept <- x$intercept[lam]
    }
  }

  if(is.matrix(X)){
    eta <- X %*% x$beta[,lam] + intercept
  }
  if(is.vector(X)){
    eta <- sum(X * x$beta[,lam]) + intercept
  }


  if(x$type == "linear"){
    y.pred <- eta
  }

  if(x$type == "logit"){
    y.pred = exp(eta)/(1+exp(eta))
  }

 if(x$type == "cox"){
    y.pred = exp(eta)
  }

  return(y.pred)
}
