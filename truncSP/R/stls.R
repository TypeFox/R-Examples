funcval.STLS <- function(bet,x,y) 		
{
	funcvalue <- sum((y-pmax((1/2)*y, x%*%bet))^2)
	return(funcvalue)
}

stls.fit<-function (formula, mf, point, direction, bet, ...) 
{
  y <- model.response(mf, "numeric")
  y <- matrix(y)
  x <- model.matrix(formula, mf)
  
    	dots <- list(...)
    	if (is.null(dots$control)) 
        	control <- list(maxit = 2000)
    	else control <- dots$control
    	z <- optim(par = bet, fn = funcval.STLS, control = control, 
        	x = x, y = y)
	    z$counts <- z$counts[1]
    	if (z$counts < 10) 
        	warning("Convergence reached after very few iterations. This could be because of incorrect \n specifications (point, direction, starting values etc.) in the function call.")
   	  b <- z$par
    	
    	z$residuals <- y - x %*% b
    	
    	if (direction == "right")
      { 
        	z$par <- -z$par
        	bet <- -bet
     	}
    
    	if (point != 0)
      { 
        	z$par[1, 1] <- z$par[1, 1] + point
          bet[1, 1] <- bet[1, 1] + point
      }
    
    	z$startcoef <- bet
    	rownames(z$startcoef) <- c(dimnames(x)[[2]])
    	colnames(z$startcoef) <- c("")
    	dimnames(z$par) <- dimnames(z$startcoef)
    	z$coefficients <- t(z$par)
    	z$df.residual <- length(mf[, 1]) - ncol(x)
    	z$fitted.values <- x %*% z$par
    	return(z)
}



stls <- function (formula, data, point = 0, direction = "left", beta = "ml", covar = FALSE, na.action, ...) 
{
    	cll <- match.call()
    	mf <- match.call(expand.dots = FALSE)
    	m <- match(c("formula", "data", "na.action"), names(mf), 0L)
    	mf <- mf[c(1L, m)]
    	mf$drop.unused.levels <- TRUE
    	mf[[1L]] <- as.name("model.frame")
    	mf <- eval(mf, parent.frame())
    	if (point != 0)
        	mf[, 1] <- mf[, 1] - point
    	
    	if (direction == "right") 
        	mf[, 1] <- -mf[, 1]
   	 
    	
	if (is.numeric(beta)) 
	{
	 bet <- matrix(beta)
 	  if (point != 0) 
    {
   	  bet[1, 1] <- bet[1, 1] - point
   	}
    if (direction == "right") 
    {
   	  bet <- -bet
   	}
	}
	else 
	{
		if (beta == "ols") 
		{
			bet <- lm(data = mf)$coef
			bet <- matrix(bet)
		}
		if (beta == "ml") 
		{
			mlcof <- mlcoef(mf)
			bet <- matrix(mlcof[-length(mlcof)])
  	}
  	if (beta != "ols" && beta != "ml") 
  	 stop("'beta' must be numeric or a valid method.")
	}

	z <- stls.fit(formula, mf, point, direction, bet, ...)
	
	if (covar) 
	{
		dots <- list(...)
		if (is.null(dots$R)) 
			R <- 2000
		else R <- dots$R
		if (is.null(dots$control)) control <- list(maxit=2000) else control <- dots$control
		bootobj <- covar.bootSTLS(mf, funcSTLS, R = R, formula, beta, bet, point, direction, control)
		z$covariance <- bootobj[[1]]
		rownames(z$covariance) <- rownames(z$startcoef)
		colnames(z$covariance) <- rownames(z$startcoef)
		z$bootrepl <- bootobj[[2]]
		z$R <- R
	}
	class(z) <- c("stls")
	z$call <- cll
z
}
