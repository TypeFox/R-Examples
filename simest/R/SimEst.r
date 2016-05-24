sim.est <- function(x, y, method = c("cvx.pen", "cvx.lip", 
                        "cvx.lse", "smooth.pen"), lambda = NULL, 
                        beta.init = NULL, w = NULL, nmulti = 5, 
                        bin.tol = 1e-04, beta.tol = 1e-04, beta.iter = 100, ...) UseMethod("sim.est")

sim.est.default <- function(x, y, method = c("cvx.pen", "cvx.lip", 
                        "cvx.lse", "smooth.pen"), lambda = NULL, 
                        beta.init = NULL, w = NULL, nmulti = 5, 
                        bin.tol = 1e-04, beta.tol = 1e-04, beta.iter = 100, ...){
	if(length(method) > 1L)
		stop("'method' should be a vector of length 1!")
	if(!(method %in% c("cvx.pen", "cvx.lip", "cvx.lse", "smooth.pen")))
		stop("The input for argument 'method' is unrecognized!")
	if(length(lambda) >= 1L && !any(method == c("cvx.pen", "cvx.lip", "smooth.pen"))){
		warning("Tuning parameter 'lambda' can only be used with 'cvx.pen', 'cvx.lip' or 'smooth.pen'!")
		method <- 'cvx.pen'
	}
	if(method == "cvx.lip" & is.null(lambda)){
		stop("The input method 'cvx.lip' requires a \ntuning parameter specification.")
	}
	x <- as.matrix(x)
	y <- as.vector(y)	
    if(nrow(x) != length(y))
        stop("Number of rows of 'x' must match the length of 'y'")	
	if (!all(is.finite(cbind(x, y)))) 
        stop("missing or infinite values in inputs are not allowed")
	n <- length(y)
	if(n <= 2)
        stop("Number of samples must be greater than 2.")
    w <- if (is.null(w)) 
        rep_len(1, n)
      else {
        if (n != length(w)) 
            stop("lengths of 'x' and 'w' must match")
        if (any(w < 0)) 
            stop("all weights should be non-negative")
        if (all(w == 0)) 
            stop("some weights should be positive")
        (w * sum(w > 0))
      }
    if(!is.null(beta.init))
        nmulti <- 1  
    if(is.null(lambda) & method != "cvx.lip")
        lambda <- 0.01*n^{-0.4}
    regress <- function(t, z, q){
    	if(method == "cvx.pen")
    		return(cvx.pen.reg(t, z, lambda, w = q))
    	if(method == "cvx.lip")
    		return(cvx.lip.reg(t, z, w = q, L = lambda))
    	if(method == "cvx.lse")
    		return(cvx.lse.reg(t, z, w = q))
    	if(method == "smooth.pen")
    		return(smooth.pen.reg(t, z, lambda, w = q))		
    }
    norm <- function(tt){sqrt(sum(tt*tt))}  
    beta.init <- if (is.null(beta.init)){
    			beta.init <- runif(ncol(x), -1, 1)
    			beta.init <- beta.init/norm(beta.init)
    			beta.init <- beta.init*{2*{beta.init[1] > 0} - 1}
    		} else{
    			if(length(beta.init) != ncol(x))
    				stop("length of 'beta.init' should match the number of columns of 'x'!")
    			if(beta.init[1] < 0 || !isTRUE(all.equal(norm(beta.init), 1)))
    				stop("'beta.init' should have norm 1 and non-negative first coordinate!")
    			(beta.init*{2*{beta.init[1] > 0} - 1}/norm(beta.init))		
    		}
    b0 <- beta.init		
    BetaPath <- rep(list(matrix(NA,nrow = beta.iter,ncol = ncol(x))),nmulti)
    ObjValPath <- rep(list(rep(NA,beta.iter)),nmulti)
    flag <- rep(0,nmulti)
    beta.conv <- rep(1,nmulti)
    iter <- rep(1,nmulti)
    Min <- rep(0,nmulti)
    for(k in 1:nmulti){
      if(k > 1){
          cat("\r multistart ",k-1," of ",nmulti," done!")
          beta.init <- runif(ncol(x), -1, 1)
          beta.init <- beta.init/norm(beta.init)
          beta.init <- beta.init*{2*{beta.init[1] > 0} - 1}
      }
      i <- 1
      while(flag[k] == 0 & i <= beta.iter){
    	t <- x%*%beta.init
    	DataMat <- cbind(t, y, x)
    	tmp <- fastmerge(DataMat, w = w, tol = bin.tol)
    	DataMat <- tmp$DataMat
    	q <- tmp$w
    	t <- DataMat[,1]
    	z <- DataMat[,2]
    	tmp <- regress(t, z, q)
    	BetaPath[[k]][i,] <- unname(beta.init)
        ObjValPath[[k]][i] <- unname(tmp$minvalue)
        deriv <- tmp$deriv
        res <- tmp$residuals
        x.mat <- DataMat[,-c(1,2)]
        b1 <- beta.init + 2*colMeans(res*deriv*x.mat)
    	b1 <- b1*{2*{b1[1] > 0} - 1}/norm(b1)
    	if(norm(b1 - beta.init) < beta.tol){
    		flag[k] <- 1
    		beta.conv[k] <- 0
    		iter[k] <- i
    	}
    	beta.init <- b1
        i <- i + 1
      }
      BetaPath[[k]] <- as.matrix(unname(as.data.frame(na.omit(BetaPath[[k]]))))
      ObjValPath[[k]] <- as.numeric(na.omit(ObjValPath[[k]]))
      if(flag[k] == 0){ iter[k] <- beta.iter; beta.conv[k] <- 1 }
      ttmp <- ObjValPath[[k]]
      Min[k] <- min(ttmp)
    }
    cat("\rmultistart ",nmulti," of ",nmulti," done!\r\n")
    K <- which.min(Min)
    B <- BetaPath[[K]][which.min(ObjValPath[[K]]),]
    t <- x%*%B
    DataMat <- unname(cbind(t, y, x))
    # print(cbind(DataMat,w))
    tmp <- fastmerge(DataMat, w = w, tol = bin.tol)
    q <- tmp$w
    # print(cbind(tmp$DataMat,q))
    t <- tmp$DataMat[,1]
    z <- tmp$DataMat[,2]
    out <- regress(t, z, q)
    ret <- list(beta = B, beta.init = b0, nmulti = nmulti, 
    	minvalue = out$minvalue, x.values = out$x.values, y.values = out$y.values)
    ret$fit.values <- out$fit.values
    ret$x.mat <- tmp$DataMat[,-c(1,2)]
    ret$BetaPath <- BetaPath
    ret$ObjValPath <- ObjValPath
    ret$convergence <- c(beta.iter = round(mean(beta.conv)), func.iter = out$convergence)
    ret$deriv <- out$deriv
    ret$residuals <- out$residuals
    ret$itervec <- iter
    ret$iter <- sum(iter)
    ret$regress <- out
    ret$method <- method
    ret$call <- match.call()
    class(ret) <- "sim.est"
    return(ret)
}

print.sim.est <- function(x,...){
	cat("Call:\n")
	print(x$call)
	cat("Initial Value of beta is:\n")
	print(x$beta.init)
	cat("Estimate of beta is:\n")
	print(x$beta)
	cat("Minimum Criterion Value Obtained:\n")
	print(x$minvalue)
	cat("Number of Iterations:\n")
	print(x$iter)
	cat("Convergence Status:\n")
	print(x$convergence)
}

plot.sim.est <- function(x,...){
	xx <- x$x.values
    yx <- x$y.values
    fitx <- x$fit.values
    resx <- x$residuals
    diagnostics = TRUE
    if(x$method == "cvx.pen")
    	tt <- "Convex Link Function Estimate\n using Penalized LS Objective"
    if(x$method == "cvx.lse")
    	tt <- "Convex Link Function Estimate\n using Least Squares Objective"
    if(x$method == "cvx.lip")
    	tt <- "Convex Link Function Estimate\n using Lipschitz LS Objective"
    if(x$method == "smooth.pen")
    	tt <- "Unconstrained Link Function Estimate\n using Penalized LS Objective"	
    if(diagnostics){
    	plot.window(c(0,7), c(0,7))
    	par(mfrow = c(2,2), mar = c(3,3,3,1), mgp = c(1.8,0.5,0))
    	plot(xx, yx, xlab = expression(x^T*hat(beta)), ylab = expression(paste('y and ',hat(y),' values')),
    		type = 'p', pch = 20, cex = 1, main = tt)
    	lines(xx, fitx, lwd = 2)
    	plot(fitx,resx,xlab = 'Fitted Values',ylab = "Residuals",pch = 20, type = 'p', main = "Fitted vs Residuals")
	    abline(h = 0.0, lty = 4)
        ttt <- max(unlist(lapply(x$ObjValPath,FUN=max)))
        tttp <- min(unlist(lapply(x$ObjValPath,FUN=min)))
        tttt <- max(x$itervec)
	    plot(seq_len(x$itervec[1]),x$ObjValPath[[1]][seq_len(x$itervec[1])],xlab = 'index',ylab = 'Objective Function Value',
	    	main = "Progression of (multistart)\n Iterative Algorithm", type = 'l', lwd = 2, ylim = c(tttp,ttt), xlim = c(1,tttt))
        if(x$nmulti > 1){
            for(i in 1:min(4,{x$nmulti-1})){
                lines(seq_len(x$itervec[i+1]),x$ObjValPath[[i+1]][seq_len(x$itervec[i+1])], lty = i)
            }
        }    
	    qqnorm(resx,main="Normal Q-Q Plot: Residuals",pch = 20)
	    qqline(resx)
    } else{
    	plot.window(c(0,7), c(0,7))
    	par(mfrow = c(1,1), mar = c(3,3,3,1), mgp = c(1.3,0.5,0))
    	plot(xx, yx, xlab = 'x',ylab = expression(paste('y and ',hat(y),' values')),
    		type = 'p', pch = "*", cex = 1, main = "Regress")
    }
}

predict.sim.est <- function(object, newdata = NULL, ...){
	req <- object$regress
	B <- object$beta
	if(!is.numeric(newdata)){
		stop("'newdata' should be numeric!")
	}
	if(!is.null(newdata)){
		newdata <- c(newdata)
		newdata <- matrix(newdata, ncol = length(B))
		t <- newdata%*%B
		return(predict(req, t))
	} else{
		warning("No 'newdata' found and so using input 'x' values")
		return(object$fit.values)
	}
}