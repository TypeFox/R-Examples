meansq <- function(a) sum(a^2)/length(a)

findnames <- function(cnames, path) {
	len <- length(path)
	name <- NULL
	for (i in 1:len) {
		name <- c(name, cnames[abs(path[i])])
	}	
}

weight.calc <- function(x, y, method) {
	n <- dim(x)[1]
	p <- dim(x)[2]
	if (method == "nsealasso") {
		lar <- lars(x, y, normalize = FALSE)
		lasso.preselect <- preselect(lar, n, p)
		beta.ols <- solve(t(x)%*%x)%*%t(x)%*%y
		sigma <- sqrt(sum((y-x%*%beta.ols)^2)/(n-p))  ### sigma can be any constant here
		se <- sqrt(diag(solve(t(x)%*%x)*sigma^2))
		or <- order(se, decreasing = TRUE)
		weight <- new.weight(lasso.preselect, se, or, beta.ols, p)
		
	} else if (method == "sealasso") {
		beta.ols <- solve(t(x)%*%x)%*%t(x)%*%y
		sigma <- sqrt(sum((y-x%*%beta.ols)^2)/(n-p))
		se <- sqrt(diag(solve(t(x)%*%x)*sigma^2))
		weight <- as.numeric(se / abs(beta.ols))
		
	} else if (method == "olsalasso") {
		beta.ols <- solve(t(x)%*%x)%*%t(x)%*%y
		weight <- as.numeric(1 / abs(beta.ols))
		
	} else if (method == "lasso") {
		weight <- 	rep(1, dim(x)[2])
		
	} 
	return(weight)
}

preselect <- function(lar, n, p) {
	RSStb <- summary(lar)
	path <- as.numeric(lar$actions)
	RSStb$BIC <- log(RSStb[,2]/n)+log(n)/n*RSStb[,1]
	step <- which(RSStb$BIC[-1]==min(RSStb$BIC[-1]))
	preselected <- rep(0,p)
	for (i in 1:step) {
		if (path[i]>0) preselected[path[i]] <- 1
		else preselected[-path[i]] <- 0
	}
	return(preselected)	
}

### assign the weight for NSEA-lasso based on the preliminary model
new.weight <- function(preselected, se, or, beta, p) {
	beta <- as.numeric(beta)
	W <- rep(0,p)
	se.front <- 1
	se.end <- p
	i <- 0
	for (i in 1:p) {
		eva <- or[i]
		if(preselected[eva]==1) {
			se.pick <- or[se.end] 
			se.end <- se.end - 1
		} else {
			se.pick <- or[se.front]
			se.front <- se.front + 1
		}
		W[eva] <-  se[se.pick] / beta[eva]
	}
	return(abs(W))
}

sealasso <- function(x, y, method = c("nsealasso", "sealasso", "olsalasso", "lasso")) {
	call <- match.call()
	method <- match.arg(method)
	TYPE <- switch(method,
					 nsealasso = "NSEA-lasso",
					 sealasso = "SEA-lasso",
					 olsalasso = "OLS-adaptive lasso",
					 lasso = "Lasso")
	
	n <- dim(x)[1]
	p <- dim(x)[2]
	cnames <- colnames(x)
	if (p < 2) stop("The number of candidate predictors must be greater than 1.")
	if (n < p) stop("Sample size must be greater than the number of predictors.")
	
					 
	### Standardize the design matrix and response
	I <- matrix(0,nrow=n,ncol=n)
	I[row(I)==col(I)] <- 1
	one <- rep(1,n)
	H <- I-one%*%t(one)/n
	D <- solve(diag(sqrt(apply(x,2,meansq))))
	x <- H %*% x %*% D
	y <- y - mean(y)
	
	### Find the condition index
	XTX <- t(x)%*%x
	eigenvalue <- eigen(XTX)$values
	eigenmin <- min(eigenvalue)
	eigenmax <- max(eigenvalue)
	index <- round(log(eigenmax/eigenmin), 1)
	
	### Find the solution path and BIC value
	weight <- weight.calc(x, y, method)
	x.weight <- x%*%solve(diag(weight))
	object <- lars(x.weight, y, normalize = FALSE)
	path <- as.numeric(object$actions)
	RSStb <- summary(object)
	BIC <- log(RSStb[,2]/n)+log(n)/n*(RSStb[,1]-1)
	Df <- RSStb[,1]
	mat <- cbind(path, Df[-1]-1, BIC[-1])

	if (length(cnames) > 0) {
		stepnames <- findnames(cnames, path)
		rownames(mat) <- stepnames
	} else {
		rownames(mat) <- 1:length(path)	
	}
	colnames(mat) <- c("Path", "Df", "BIC")
	
	### Find the corresponding estimated coefficients transformed back to the original scale
	## optim.step <- which.min(BIC[-1])
	beta <- object$beta[-1,]
	beta <- beta%*%D
	colnames(beta) <- cnames
	rownames(beta) <- 1:length(path)
	
	### Find the optimal estimated coefficients
	optim.step <- which.min(BIC[-1])
	optim.beta <- object$beta[optim.step+1,]%*%D
	colnames(optim.beta) <- cnames
	rownames(optim.beta) <- optim.step
	
	lst <- list(call = call, method = TYPE, weight = weight, condition.index = index, path = mat, beta = beta, optim.beta = optim.beta)
	class(lst) <- "sealasso"
	return(lst)			
}

summary.sealasso<-function(object, ...){  
  lst <- list(method = object$method, condition.index = object$condition.index, optim.beta = object$optim.beta)
  return (lst)
}


