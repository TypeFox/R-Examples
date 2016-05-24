cvx.pen.reg <- function(t, z, lambda, w = NULL, max.iter = 500, alpha.tol= 1e-03, ...) UseMethod("cvx.pen.reg")

cvx.pen.reg.default <- function(t, z, lambda, w = NULL, max.iter = 500, alpha.tol= 1e-03, ...){
  if(length(t) != length(z))
    stop("'x' and 'y' must have same length.")  
  if (!all(is.finite(c(t, z)))) 
    stop("missing or infinite values in inputs are not allowed")
  n <- length(t)
  if(n <= 2)
    stop("Number of samples must be greater than 2.")
  w <- if (is.null(w)) 
        rep_len(1, n)/n
    else {
        if (n != length(w)) 
            stop("lengths of 'x' and 'w' must match")
        if (any(w <= 0)) 
            stop("all weights should be positive")
        w/sum(w > 0)
    }
  q <- w
  p <- lambda
  tol <- alpha.tol
  if (!is.finite(p) || p < 0)
    stop("'lambda' must be non-negative and finite")  
  A <- cbind(t, z, w)
  A <- A[order(A[,1]),]
  t <- A[,1]
  z <- A[,2]
  w <- A[,3]
  n <- length(t)
  N <- n-2
  K <- matrix(0,nrow = N, ncol = N+2)
  a.i <- diff(t,2)/4
  for(i in 1:N){
    K[i,i] <- 1/{{t[i+1]-t[i]}*{t[i+2]-t[i]}}
    K[i,i+1] <- -1/{{t[i+2]-t[i+1]}*{t[i+1]-t[i]}}
    K[i,i+2] <- 1/{{t[i+2]-t[i]}*{t[i+2]-t[i+1]}}
  }
  a.init <- a.i
  Mfun <- function(a){
    mf <- .C("Mfun", as.integer(N), as.double(a), as.double(t), double(2*N + 2), 
    	double(N*N), PACKAGE = "simest")
    Vmat <- mf[[4]]
    V <- matrix(Vmat, ncol = 2, nrow = N+1, byrow = TRUE)
    Mmat <- mf[[5]]
    M <- matrix(Mmat, ncol = N, nrow = N, byrow = TRUE)
    return(list(M = M, V = V))
  }
  norm = function(a){
    return(sqrt(sum(a^2)))
  }
  M11 <- p*K%*%(t(K)/q) 
  M22 <- diff(diff(z)/diff(t))/diff(t,2)
  ba <- a.init
  MM <- Mfun(ba)$M
  bf <- norm({{MM + M11}%*%ba} - M22)
  run <- function(a0, sup.iter = max.iter){
    loop <- rep(0,sup.iter/50)
    for(i in 1:{sup.iter}){
      MF <- Mfun(a0)
      Bf0 <- MF$M + M11
      a1 <- solve.pentadiag(Bf0, M22)
      check <- norm(Bf0%*%a0 - M22)
      if(check < bf){
        ba <- a0
        bf <- check
      }
      if(i%%50 == 0){
        if(is.element(check,loop)){
          #print(paste('Code Looped at stop criterion value',check))
          #warning('Fixed point iteration started looping.\n And so returning the best value obtained.')
          flag <- 2
          temp <- as.vector(p*(t(K)/q)%*%ba)
          temp1 <- t(ba)%*%M22
          zhat <- z - temp
          V <- MF$V
          DD <- .C("deriv", as.integer(n), as.double(zhat), as.double(ba), 
          	as.double(c(t(V))), as.double(t), double(n), PACKAGE = "simest")
          D <- DD[[length(DD)]]
          return(list(alpha = ba, minvalue = c(p*temp1), x.values = t, y.values = z, 
            fit.values = zhat, residuals = temp, convergence = flag, 
            Kz = M22, Vmat = MF$V, deriv = D, iter = i))
        }
        else loop[i/50] <- check
      }
      if(check <= tol){
        temp <- as.vector(p*(t(K)/q)%*%a1)
        temp1 <- t(a1)%*%M22
        zhat <- z - temp
        flag <- 0
        V <- Mfun(a1)$V
        DD <- .C("deriv", as.integer(n), as.double(zhat), as.double(a1),
        	as.double(c(t(V))), as.double(t), double(n), PACKAGE = "simest")
        D <- DD[[length(DD)]]
        return(list(alpha = a1, minvalue = c(p*temp1), x.values = t, y.values = z, 
          fit.values = zhat, residuals = temp, convergence = flag,
          Kz = M22, Vmat = MF$V, deriv = D, iter = i))
      }
      if(i == sup.iter){
        #warning('sup.iter in convex function estimation reached.')
        flag <- 1
        temp <- as.vector(p*(t(K)/q)%*%ba)
        temp1 <- t(ba)%*%M22
        zhat <- z - temp
        V <- MF$V
        DD <- .C("deriv", as.integer(n), as.double(zhat), as.double(ba),
        	as.double(c(t(V))), as.double(t), double(n), PACKAGE = "simest")
        D <- DD[[length(DD)]]
        return(list(alpha = a1, minvalue = c(p*temp1), x.values = t, y.values = z, 
          fit.values = zhat, residuals = temp, convergence = flag, 
          Kz = M22, Vmat = MF$V, deriv = D, iter = i))
      }
      else a0 <- a1
      i <- i + 1
    }
  }
  ret <- run(a.init)
  ret$call <- match.call()
  class(ret) <- "cvx.pen.reg"
  return(ret)
}

print.cvx.pen.reg <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("Minimum Criterion Value Obtained:\n")
  print(x$minvalue)
  cat("Number of Iterations:\n")
  print(x$iter)
  cat("Convergence Status:\n")
  print(x$convergence)  
  #cat("Input 'x', 'y', the fit and the derivative:\n")
  #write.table(format(rbind(c('x','y','yhat','deriv'),
  #  cbind(object$x.values, object$y.values, object$fit.values, object$deriv)), 
  #  justify="right"), row.names=F, col.names = F, quote=F)
}

plot.cvx.pen.reg <- function(x, ...){
  xx <- x$x.values
  yx <- x$y.values
  fitx <- x$fit.values
  resx <- x$residuals
  diagnostics = TRUE
  if(diagnostics){
    plot.window(c(0,7), c(0,7))
    par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(1.3,.5,0))
    plot(xx,yx,xlab = 'x',ylab = expression(paste('y and ',hat(y),' values')), 
      type = 'p', pch = "*", cex = 1, main = "Convex Regression using\n Penalized Least Squares")
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
      type = 'p', pch = "*", cex = 1, main = "Convex Regression using\n Penalized Least Squares")
    lines(xx, fitx, lwd = 2)    
  }
  invisible(list(x = xx, y = yx, fit = fitx))
}

predict.cvx.pen.reg <- function(object, newdata = NULL, ...){
  t <- unname(object$x.values)
  zhat <- unname(object$fit.values)
  a <- unname(object$alpha)
  D <- unname(object$deriv)
  V <- unname(object$Vmat)
  V <- c(t(V))
  n <- length(t)
  if(is.null(newdata)){
    warning("No 'newdata' found and so using input 'x' values")
    return(zhat)
  } else{
    newdata <- as.vector(newdata)
    r <- length(newdata)
    dim <- c(n,r)
    out <- .C("pred", as.integer(dim), as.double(t), as.double(zhat), as.double(a),
      as.double(D), as.double(V), as.double(newdata), PACKAGE = "simest")
    return(out[[7]])
  }  
}