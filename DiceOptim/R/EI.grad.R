EI.grad <- function(x, model, plugin=NULL, type="UK", minimization = TRUE, envir=NULL){ 
	
  ########################################################################################
  if (is.null(plugin)){ 
    if (minimization) {
      plugin <- min(model@y)
    } else {
      plugin <- -max(model@y)
    }
  }
  m <- plugin
  
  ########################################################################################
  # Convert x in proper format(s)
  d <- length(x)
  if (d != model@d){ stop("x does not have the right size") }
  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) = colnames(model@X)
  ########################################################################################
  # Get quantities related to the model
  T <- model@T
  X <- model@X
  z <- model@z
  u <- model@M
  covStruct <- model@covariance
  
  # Get quantities related to the prediction
  if (is.null(envir))
  {  
    predx <- predict(object=model, newdata=newdata, type=type, checkNames = FALSE,se.compute=TRUE,cov.compute=FALSE)
     kriging.mean <- predx$mean
     if(!minimization) kriging.mean <- -kriging.mean
     kriging.sd <- predx$sd
     v <- predx$Tinv.c
     c <- predx$c
  
     xcr <- (m - kriging.mean)/kriging.sd
     xcr.prob <- pnorm(xcr)
     xcr.dens <- dnorm(xcr)    
  } else
  {  # If uploaded through "envir", no prediction computation is necessary 
     toget <- matrix(c("xcr", "xcr.prob", "xcr.dens", "kriging.sd", "c", "Tinv.c"), 1, 6)
     apply(toget, 2, get, envir=envir)
     xcr        <- envir$xcr
     xcr.prob   <- envir$xcr.prob
     xcr.dens   <- envir$xcr.dens
     kriging.sd <- envir$kriging.sd
     c          <- envir$c
     v          <- envir$Tinv.c
  }

  F.newdata <- model.matrix(model@trend.formula, data=newdata)

  ########################################################################################
  # Pursue calculation only if standard deviation is non-zero
  if ( kriging.sd/sqrt(model@covariance@sd2) < 1e-06) 
  { ei.grad <- rep(0,d)
  } else 
  { # Compute derivatives of the covariance and trend functions
    dc <- covVector.dx(x=newdata.num, X=X, object=covStruct, c=c)  
    f.deltax <- trend.deltax(x=newdata.num, model=model)
  
    # Compute gradients of the kriging mean and variance
    W <- backsolve(t(T), dc, upper.tri=FALSE)
    kriging.mean.grad <- t(W)%*%z + t(model@trend.coef%*%f.deltax)
    if (!minimization) kriging.mean.grad <- -kriging.mean.grad
    if (type=="UK")
    { tuuinv <- solve(t(u)%*%u)
      kriging.sd2.grad <-  t( -2*t(v)%*%W +
                            2*(F.newdata - t(v)%*%u )%*% tuuinv %*%
                            (f.deltax - t(t(W)%*%u) ))
    } else
    { kriging.sd2.grad <-  t( -2*t(v)%*%W) }
    
    kriging.sd.grad <- kriging.sd2.grad / (2*kriging.sd)

    # Compute gradient of EI
    ei.grad <- - kriging.mean.grad * xcr.prob + kriging.sd.grad * xcr.dens
  }
  ########################################################################################  
  return(ei.grad)
}