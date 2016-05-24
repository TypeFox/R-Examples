bothsidesmodel.mle <-
function(x,y,z=diag(qq),pattern=matrix(1,nrow=p,ncol=l)) {
	x <- cbind(x)
	y <- cbind(y)
	n <- nrow(y)
	p <- ncol(x)
	qq <- ncol(y)
	z <- cbind(z)
	l <- ncol(z)

	if(length(x)==1&&x==0) {
		bsm <- list(Beta=0,SE=0,T=0,Covbeta=0,df=n)
		resid <- y
		psum <- 0
	} else {
		bsm <- list()
		beta <- se <- tt <- matrix(0,nrow=p,ncol=l)
		psum <- sum(pattern)
		if(psum>0) {
			rowin <- apply(pattern,1,max)==1
			colin <- apply(pattern,2,max)==1
			x1 <- x[,rowin,drop=FALSE]
			z1 <- z[,colin,drop=FALSE]
			yz <- cbind(y%*%solve(t(fillout(z1))))
			l1 <- ncol(z1)
			p1 <- ncol(x1)
			yza <- yz[,1:l1,drop=FALSE]
			xyzb <- x[,rowin,drop=FALSE]
			if(l1<qq) xyzb <- cbind(xyzb,yz[,(l1+1):qq])
			pattern1 <- pattern[rowin,colin]
			if(min(pattern1)==1) {
				bsm <- bsm.simple(xyzb,yza,diag(l1))
				bsm$Cx <- bsm$Cx[1:p1,1:p1]
			} else {
				pattern2 <- rbind(pattern1,matrix(1,nrow=qq-l1,ncol=l1))
				bsm <- bsm.fit(xyzb,yza,diag(l1),pattern2)
			}
			beta[rowin,colin] <- bsm$Beta[1:p1,]
			se[rowin,colin] <- bsm$SE[1:p1,]
			tt[rowin,colin] <- bsm$T[1:p1,]
			bsm$Covbeta <- bsm$Covbeta[1:psum,1:psum]
		}			
		bsm$Beta <- beta
		bsm$SE <- se
		bsm$T <- tt
		if(psum==0) {
			bsm$Covbeta <- 0
			bsm$df <- n
			p1 <- 0
		}
		resid <- y-x%*%beta%*%t(z)
	}
	bsm$SigmaR <- t(resid)%*%resid/n
	bsm$Deviance <- n*logdet(bsm$SigmaR) + n*qq
	bsm$Dim <- psum+qq*(qq+1)/2
	bsm$AICc <- bsm$Deviance+(n/(n-p1-qq-1))*2*bsm$Dim
	bsm$BIC <- bsm$Deviance + log(n)*bsm$Dim
	bsm$Cx <- NULL
	bsm$Sigmaz <- NULL
	
	bsm
}
