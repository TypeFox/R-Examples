bothsidesmodel <-
function(x,y,z=diag(qq),pattern=matrix(1,nrow=p,ncol=l)) {
	x <- cbind(x)
	y <- cbind(y)
	n <- nrow(y)
	p <- ncol(x)
	qq <- ncol(y)
	z <- cbind(z)
	l <- ncol(z)
	
	if((p*nrow(x)==1&&x==0) || (l*nrow(z)==1&&z==0)) {
		sigma <- t(y)%*%y/n
		output <- list(Beta = 0, SE = 0,T = 0, Covbeta = 0, df = n, Sigmaz=sigma)
		return(output)
	}

	yz <- t(lm(t(y)~z-1)$coef)
	residb <- lm(yz~x-1)$resid
	sigmab <- t(residb)%*%residb/(n-p)

	if(sum(pattern)==0) {
		residss <- ((n-p-l-1)/(n-p))*tr(solve(sigmab,t(yz)%*%yz))
		output <- list(Beta = 0, SE = 0,T = 0, Covbeta = 0, df = n,Sigmaz=sigmab, 
		ResidSS=residss,Dim=l*(l+1)/2,Cp=residss+l*(l+1))
		return(output)
	}

	
	rowse <- rowt <- rowb <- rep(0,p*l)
	rowp <- c(t(pattern))==1
	
	xx <- t(x)%*%x
	xyz <- t(x)%*%yz
	xyz <- c(t(xyz))[rowp]
	xxzz <- kronecker(xx,diag(l))[rowp,rowp]
	dstarinv <- solve(xxzz)
	g <- xyz%*%dstarinv
	rowb[rowp] <- g
	beta <- matrix(rowb,nrow=p,byrow=T)
	
	df <- bothsidesmodel.df(xx,n,pattern)
	residz <- yz-x%*%beta
	residsscp <- t(residz)%*%residz
	sigmaz <- residsscp/df
	xxzsz <- kronecker(xx,sigmaz)[rowp,rowp]
	covbeta <- dstarinv%*%xxzsz%*%dstarinv
	covbetab <- matrix(0,p*l,p*l)
	covbetab[rowp,rowp] <- covbeta
	se <- sqrt(diag(covbeta))
	tt <- g/se
	rowse[rowp] <- se
	rowt[rowp] <- tt	
	residss <- ((n-p-l-1)/(n-p))*tr(solve(sigmab,residsscp))
	dd <- sum(pattern)+l*(l+1)/2

	
	list(Beta = beta, SE = matrix(rowse,nrow=p,byrow=T), 
		T = matrix(rowt,nrow=p,byrow=T), Covbeta = covbetab, df = diag(df),Sigmaz = sigmaz,
		ResidSS = residss,Dim=dd,Cp=residss+2*dd)
}
