AEI.grad <- function(x, model, new.noise.var=0, y.min=NULL, type = "UK", envir=NULL){

  ########## Convert x in proper format(s) ###
  d <- length(x)
  if (d != model@d){ stop("x does not have the right size") }
  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) = colnames(model@X)
  
  # Compute y.min if missing
  if (is.null(y.min))
  {
    pred <- predict(object=model, newdata=model@X, type=type, checkNames = FALSE)
    mk <- pred$mean
    sk <- pred$sd
    qk <- mk + qnorm(0.75)*sk
    y.min <- mk[which.min(qk)]
  }
  
  T <- model@T
  X <- model@X
  z <- model@z
  u   <- model@M
  covStruct <- model@covariance
  
  if (is.null(envir))
  {  predx <- predict.km(object=model, newdata=newdata, type=type, checkNames = FALSE)
     mk  <- predx$mean
     sk  <- predx$sd
     c   <- predx$c
     v   <- predx$Tinv.c
     xcr <- (y.min - mk)/sk
     xcr.prob <- pnorm(xcr)
     xcr.dens <- dnorm(xcr)
  } else
  { toget <- matrix(c("xcr", "xcr.prob", "xcr.dens", "c", "Tinv.c", "mk","sk"),1,7)
    apply(toget, 2, get, envir=envir)
    xcr      <- envir$xcr
    xcr.prob <- envir$xcr.prob
    xcr.dens <- envir$xcr.dens
    c        <- envir$c
    v        <- envir$Tinv.c
    mk       <- envir$mk
    sk       <- envir$sk
  }
  F.newdata <- model.matrix(model@trend.formula, data=newdata)
  
  ###########################################################################
  
  if (sk < sqrt(model@covariance@sd2)/1e6)
  { aei.grad.val <- rep(0,d)
  } else 
  { # AEI
    ei.val <- (y.min - mk) * xcr.prob + sk * xcr.dens
    pen <- (1- sqrt(new.noise.var)/sqrt(new.noise.var + sk^2))
    aei.val <- ei.val * pen
    
    # EI gradient
    # Compute derivatives of the covariance and trend functions
    dc <- covVector.dx(x=newdata.num, X=X, object=covStruct, c=c)  
    f.deltax <- trend.deltax(x=newdata.num, model=model)
    
    # Compute gradients of the kriging mean and variance
    W <- backsolve(t(T), dc, upper.tri=FALSE)
    mk.grad <- t(W)%*%z + t(model@trend.coef%*%f.deltax)
    
    if (type=="UK")
    { tuuinv <- solve(t(u)%*%u)
      sk2.grad <-  t( -2*t(v)%*%W +
                      2*(F.newdata - t(v)%*%u )%*% tuuinv %*%
                      ( f.deltax - t(t(W)%*%u) ))
    } else
    { sk2.grad <-  t( -2*t(v)%*%W)
    }
    sk.grad <- sk2.grad / (2*sk)
    
    # Compute gradient of EI
    ei.grad <- - mk.grad * xcr.prob + sk.grad * xcr.dens
    
    # Penalization gradient
    pen.grad <- sqrt(new.noise.var)/2*sk2.grad*(new.noise.var + sk^2)^(-3/2)
    
    # AEI gradient
    aei.grad.val <- ei.grad*pen + ei.val*pen.grad
    
#     pen.grad <- - ei.val * (1 - sqrt(new.noise.var))*sk*sk.grad*(new.noise.var + sk^2)^(-3/2)
#     aei.grad.val <- ei.grad*pen + ei.val^2*pen.grad
  }
  return(aei.grad.val)
}
