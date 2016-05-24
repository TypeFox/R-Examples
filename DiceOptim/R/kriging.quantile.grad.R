#source("kriging.quantile.R")

kriging.quantile.grad <- function(x, model, beta=0.1, type="UK", envir=NULL)
{
  ########## Convert x in proper format(s) ###
  d <- length(x)
  if (d != model@d){ stop("x does not have the right size") }
  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) = colnames(model@X)

	T <- model@T
	X <- model@X
	z <- model@z
  u <- model@M
  
	covStruct <- model@covariance
  
	if (is.null(envir))
	{ predx <- predict(object=model, newdata=newdata, type=type, checkNames = FALSE)
	  mk <- predx$mean
	  sk <- predx$sd
	  c <- predx$c 
	  v <- predx$Tinv.c
	} else
	{ toget <- matrix(c("mk", "sk",  "c", "Tinv.c"),1,4)
	  apply(toget, 2, get, envir=envir)
	  c   <- envir$c
	  v   <- envir$Tinv.c
	  mk  <- envir$mk
	  sk  <- envir$sk
	}
  F.newdata <- model.matrix(model@trend.formula, data=newdata)
  
	# Compute derivatives of the covariance and trend functions
	dc <- covVector.dx(x=newdata.num, X=X, object=covStruct, c=c)  
	f.deltax <- trend.deltax(x=newdata.num, model=model)
	
	# Compute gradients of the kriging mean and variance
	W <- backsolve(t(T), dc, upper.tri=FALSE)
	mk.grad <- t(W)%*%z + t(model@trend.coef%*%f.deltax)
	
  if (sk < sqrt(model@covariance@sd2)/1e6)
  { sk.grad <- 0
  } else
  { 
    if (type=="UK")
    { tuuinv <- solve(t(u)%*%u)
	    sk2.grad <-  t( -2*t(v)%*%W +
	             2*(F.newdata - t(v)%*%u )%*% tuuinv %*%
	                          (f.deltax - t(t(W)%*%u) ))
    } else { sk2.grad <-  t( -2*t(v)%*%W)}
	  
	  sk.grad <- sk2.grad / (2*sk)
  }

  return(quantile.grad <- mk.grad + qnorm(beta)*sk.grad)
}
