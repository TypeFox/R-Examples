smooth.pen.reg <- function(x, y, lambda, w = NULL,...) UseMethod("smooth.pen.reg")

smooth.pen.reg.default <- function(x, y, lambda, w = NULL,...){
  if (!all(is.finite(c(x, y)))) 
    stop("missing or infinite values in inputs are not allowed")
  if(length(x) != length(y))
    stop("'x' and 'y' must have same length.")
  n <- length(x)  
  w <- if (is.null(w)) 
        rep_len(1, n)/n
    else {
        if (n != length(w)) 
            stop("lengths of 'x' and 'w' must match")
        if (any(w <= 0)) 
            stop("all weights should be positive")
        (w * sum(w > 0))/sum(w)
    }  
  A <- cbind(as.vector(unname(x)),as.vector(unname(y)), as.vector(w))
  #print(dim(A))
  A <- A[order(A[,1]),]
  x <- as.vector(A[,1])
  y <- as.vector(A[,2])
  w <- as.vector(A[,3])
  if (!is.finite(lambda) || lambda < 0)
    stop("'lambda' must be non-negative and finite")   
  h <- diff(x)
  R <- matrix(0,nrow = n-2, ncol = n-2)
  Q <- matrix(0,nrow = n, ncol = n-2)
  R[n-2,n-2] <- 2*{h[n-2] + h[n-1]}/3
  Q[n-2,n-2] <- 1/h[n-2]; Q[n-1,n-2] = -({1/h[n-2]} + {1/h[n-1]}); Q[n,n-2] = 1/h[n-1]
  for(i in 1:{n-3}){
    R[i,i] <- 2*{h[i] + h[i+1]}/3
    R[i,i+1] <- h[i+1]/3
    R[i+1,i] <- R[i,i+1]
    Q[i,i] <- 1/h[i]
    Q[i+1,i] <- -{1/h[i]} - {1/h[i+1]}
    Q[i+2,i] <- 1/h[i+1]
  }
  a <- lambda
  Qty <- diff(diff(y)/h)
  T <- R + a*t(Q)%*%(w*Q)
  C <- as.vector(solve.pentadiag(T, Qty))
  m <- as.vector(y - a*w*{Q%*%C})
  tmp <- splinefun(x, m)
  deriv <- unname(tmp(x, deriv = 1))
  ret <- list(x.values = x, y.values = y, fit.values = m, 
    deriv = as.vector(deriv), residuals = y - m, iter = 1, convergence = 0)
  ret$minvalue <- mean(w*{y-m}^2)
  ret$call <- match.call()
  ret$splinefun <- tmp
  class(ret) <- "smooth.pen.reg"
  return(ret)
}

print.smooth.pen.reg <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("Minimum Criterion Value Obtained:\n")
  print(x$minvalue)
  cat("Number of Iterations:\n")
  print(x$iter)
  cat("Convergence Status:\n")
  print(x$convergence)	
}

plot.smooth.pen.reg <- function(x, ...){
  xx <- x$x.values
  yx <- x$y.values
  fitx <- x$fit.values
  resx <- x$residuals
  diagnostics = TRUE
  if(diagnostics){
	plot.window(c(0,7), c(0,7))
	par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(1.3,.5,0))
	plot(xx,yx,xlab = 'x',ylab = expression(paste('y and ',hat(y),' values')), 
	  type = 'p', pch = "*", cex = 1, main = "Smoothing Spline using\n Penalized Least Squares")
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
	  type = 'p', pch = "*", cex = 1, main = "Smoothing Spline using\n Penalized Least Squares")
	lines(xx, fitx, lwd = 2)  	
  }
  invisible(list(x = xx, y = yx, fit = fitx))
}

predict.smooth.pen.reg <- function(object, newdata = NULL,...){
	if(!is.null(newdata))
		return(object$splinefun(as.vector(newdata)))
	else{
		warning("No 'newdata' found and so using input 'x' values")
		return(object$fit.values)
	}
}