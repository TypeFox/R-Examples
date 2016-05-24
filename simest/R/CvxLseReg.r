cvx.lse.reg <- function(t, z, w = NULL, ...) UseMethod("cvx.lse.reg")

cvx.lse.reg.default <- function(t, z, w = NULL, ...){
	t <- as.vector(t)
	z <- as.vector(z)
	if (!all(is.finite(c(t, z)))) 
    stop("missing or infinite values in inputs are not allowed")
    if(length(t) != length(z))
      stop("'x' and 'y' must have same length.")
    n = length(t)
    if(n <= 2)
      stop("Number of samples must be greater than 2.")
    w <- if (is.null(w)) 
        rep_len(1, n)
      else {
        if (n != length(w)) 
            stop("lengths of 'x' and 'w' must match")
        if (any(w < 0)) 
            stop("all weights should be non-negative")
        (w * sum(w > 0))
      }
  A <- cbind(t, z, w)  
	A <- A[order(A[,1]),]
	x <- as.vector(A[,1])
	y <- as.vector(A[,2])
  w <- as.vector(A[,3])
	A <- matrix(0,nrow = n-2,ncol = n)
	for(i in 1:{n-2}){
		A[i,i] <- x[i+2] - x[i+1]
		A[i,i+1] <- x[i] - x[i+2]
		A[i,i+2] <- x[i+1] - x[i]
	}
	G <- t(t(A)/sqrt(w))
	h <- - A%*%y
	E <- rbind(t(G), t(h))
	f <- c(rep(0,n),1)
	tmp <- nnls(E,f)
	u <- tmp$x
	r <- as.vector(E%*%u - f)
	z <- -r[1:n]/r[n+1]
	fit <- y + z/sqrt(w)
	deriv <- diff(fit)/diff(x)
	deriv <- c(deriv, deriv[length(deriv)])
	ret1 <- list(x.values = x, y.values = y, fit.values = fit, deriv = deriv, iter = 1,
		residuals = y - fit, minvalue = mean({w*{y-fit}^2}), convergence = tmp$mode)
	ret1$call <- match.call()
	class(ret1) <- "cvx.lse.reg" 
	return(ret1)
}

print.cvx.lse.reg <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("Minimum Criterion Value Obtained:\n")
  print(x$minvalue)
  cat("Number of Iterations:\n")
  print(x$iter)
  cat("Convergence Status:\n")
  print(x$convergence)
}

plot.cvx.lse.reg <- function(x, ...){
  xx <- x$x.values
  yx <- x$y.values
  fitx <- x$fit.values
  resx <- x$residuals
  diagnostics = TRUE
  if(diagnostics){
    plot.window(c(0,7), c(0,7))
    par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(1.3,.5,0))
    plot(xx,yx,xlab = 'x',ylab = expression(paste('y and ',hat(y),' values')), 
      type = 'p', pch = "*", cex = 1, main = "Convex Regression using\n Least Squares")
    lines(xx, fitx, lwd = 2)
    plot(fitx,resx,xlab = 'Fitted Values',ylab = "Residuals",pch = "*", type = 'p', main = "Fitted vs Residuals")
    abline(h = 0.0, lty = 4)
    plot(yx,fitx,xlab = "Actual Values",ylab = "Fitted Values",pch = "*", type = 'p', main = "Actual vs Fitted")
    abline(a = 0, b = 1, lty = 4)
    qqnorm(resx)
    qqline(resx)
  } else{
    plot.window(c(0,7), c(0,7))
    par(mfrow=c(1,1), mar=c(3,3,3,1), mgp=c(1.3,.5,0))
    plot(xx,yx,xlab = 'x',ylab = expression(paste('y and ',hat(y),' values')), 
      type = 'p', pch = "*", cex = 1, main = "Convex Regression using\n Least Squares")
    lines(xx, fitx, lwd = 2)    
  }
  invisible(list(x = xx, y = yx, fit = fitx))
}

predict.cvx.lse.reg <- function(object, newdata = NULL, ...){
	x <- object$x.values
	y <- object$y.values
	fit <- object$fit.values
	deriv <- object$deriv
  n <- length(x)
	if(is.null(newdata)){
    	warning("No 'newdata' found and so using input 'x' values")
    	return(fit)
  	} else{
    	newdata <- as.vector(newdata)
    	r <- length(newdata)
    	foo <- function(t){
    		if(t < x[1]){
    			return(fit[1] + deriv[1]*{t - x[1]})
    		}
    		if(t > x[n]){
    			return(fit[n] + deriv[n]*{t - x[n]})
    		}
    		for(i in 1:{n-1}){
      			if(t >= x[i] & t <= x[i+1]){
        			return(fit[i] + deriv[i]*{t - x[i]})
      			}
    		}
  		}
  		return(as.vector(unlist(sapply(newdata[seq_len(r)],foo,simplify = TRUE,USE.NAMES = FALSE))))
	}
}	