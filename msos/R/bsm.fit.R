bsm.fit <-
function(x,y,z,pattern) {
	n <- nrow(y)
	p <- ncol(x)
	l <- ncol(z)
	xx <- t(x)%*%x
	xy <- t(x)%*%y
	rowp <- c(t(pattern))==1
	lsreg <- lm(y~x-1)
	yqxy <- t(lsreg$resid)%*%lsreg$resid
	beta0 <- lsreg$coef
	sigma <- cov(y)*(n-1)/n
	rowse <- rowt <- rowb <- rep(0,p*l)
	dev0 <- n*logdet(sigma)

	maxiter = 25
	iter = 0;

	repeat {
		sigmainvz <- solve(sigma,z)
		xyz <- c(t(xy%*%sigmainvz))[rowp]
		xxzz <- kronecker(xx,t(z)%*%sigmainvz)[rowp,rowp]
		dstarinv <- solve(xxzz)
		gamma <- xyz%*%dstarinv
		rowb[rowp] <- gamma
		beta <- matrix(rowb,nrow=p,byrow=T)	
		betadiff <- beta0-beta%*%t(z)
		sigma <- (yqxy+t(betadiff)%*%xx%*%betadiff)/n
		dev <- n*logdet(sigma)	
		iter <- iter+1
		if(abs(dev-dev0) < 0.00001) break;
		if(iter>=maxiter) break;
		dev0 <- dev
	}	
	df <- n-p
	covbeta <- solve(xxzz)*(n/df)
	se <- sqrt(diag(covbeta))
	tt <- gamma/se
	rowse[rowp] <- se
	rowt[rowp] <- tt	
	se <- matrix(rowse,nrow=p,byrow=T)		
	tt <- matrix(rowt,nrow=p,byrow=T)		
	list(Beta = beta, SE = se, T = tt, Covbeta = covbeta, df = df)
}
