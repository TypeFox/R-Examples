"qrisk" <-
function(x, alpha=c(.1,.3), w = c(.7,.3), mu = .07, R = NULL, r = NULL, lambda = 10000){
#
# find optimal Choquet-risk portfolios given:
#
#	x 	(n by p) matrix of asset returns
#	alphas	alphas defining a Choquet capacity risk function
#	w	w defining weights for Choquet capacity risk function
#	R	Matrix defining constraints on the parameters 
#	r	rhs defining constraints on the parameters 
#	mu	required mean rate of return	
#	lambda	Lagrange multiplier for RoR constraint 
#
	n <- nrow(x)
	p <- ncol(x)
	m <- length(alpha)
	if(length(w)!=m)stop("length of w doesn't match length of alpha")
	xbar <- apply(x,2,mean)
	y <- x[,1]
	r <- c(r,lambda*(xbar[1]-mu), -lambda*(xbar[1]-mu))
	X <- x[,1]-x[,-1]
	R <- rbind(R,lambda*(xbar[1]-xbar[-1]), -lambda*(xbar[1]-xbar[-1]))
	R <- cbind(matrix(0,nrow(R),m),R)
	f <- rq.fit.hogg(X,y,taus=alpha,weights=w,R=R,r=r)
	fit <- f$coefficients
	pihat <- c(1-sum(fit[-(1:m)]),fit[-(1:m)])
	x <- as.matrix(x)
	yhat <- x%*%pihat
	etahat <- quantile(yhat,alpha)
	muhat <- mean(yhat)
	qrisk <- 0
	for(i in 1:length(alpha))
		qrisk <- qrisk + w[i]*sum(yhat[yhat<etahat[i]])/(n*alpha[i]) 
	list(pihat = pihat, muhat = muhat, qrisk = qrisk)
}
"srisk" <-
function(x,mu=.07,lambda=100000000,alpha=.1,eps=.0001){
#
# find optimal sigma-risk (Markowitz) portfolios given:
#
#	x 	(n by p) matrix of asset returns
#	mu	required mean return
#	lambda	Lagrange multiplier for required mean return constraint
#
	n <- nrow(x)
	p <- ncol(x)
	X <- rbind(x,apply(x,2,mean))
	y <- X[,1]
	y[n+1] <- lambda*(y[n+1]-mu)
	X <- X[,1]-X[,-1]
	X <- cbind(c(rep(1,n),0),X)
	X[n+1,] <- lambda*X[n+1,]
	fit <- lm(y~X-1)
	pihat <- c(1-sum(fit$coef[-1]),fit$coef[-1])
	if(abs(fit$residual[n+1]) > eps) return("lambda too small?")
	yhat <- x%*%pihat
	muhat <- mean(x%*%pihat)
	sigma <- sqrt(var(x%*%pihat))
	list(pihat = pihat, muhat = muhat, sigma = sigma)
}
