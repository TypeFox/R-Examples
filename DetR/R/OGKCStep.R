OGKCStep<-function(x0,scale_est,alpha=0.5){
	p<-ncol(x0)
	n<-nrow(x0)
	quanf<-function(n,p,alpha) return(floor(2*floor((n+p+1)/2)-n+2*(n-floor((n+p+1)/2))*alpha))
	initset<-function(data,scalefn,P,h){
		stopifnot(length(d <- dim(data)) == 2, length(h) == 1, h >= 1)
		n <- d[1]
		stopifnot(h <= n)
		lambda <- doScale(data %*% P, center=median, scale=scalefn)$scale
		sqrtcov    <- P %*% (lambda * t(P)) ## == P %*% diag(lambda) %*% t(P)
		sqrtinvcov <- P %*% (t(P) / lambda) ## == P %*% diag(1/lambda) %*% t(P)
		estloc <- colMedians(data %*% sqrtinvcov) %*% sqrtcov
		centeredx <- (data - rep(estloc, each=n)) %*% P
		sort.list(mahalanobisD(centeredx, FALSE, lambda))[1:h]# , partial = 1:h
	}
	ogkscatter<-function(Y,scalefn,only.P=TRUE){
		stopifnot(length(p <- ncol(Y)) == 1, p >= 1)
		U <- diag(p)

		for(i in seq_len(p)[-1L]) {# i = 2:p
		    sYi <- Y[,i]
		    ii <- seq_len(i - 1L)
		    for(j in ii) {
		        sYj <- Y[,j]
		        U[i,j] <- (scalefn(sYi + sYj)^2 - scalefn(sYi - sYj)^2) / 4
		    }
		    ## also set the upper triangle
		    U[ii,i] <- U[i,ii]
		}

		## now done above: U <- lower.tri(U) * U + t(U)    #    U <- tril(U, -1) + t(U)
		P <- eigen(U, symmetric=TRUE)$vectors
		if(only.P)
		    return(P)
		## else :
		Z <- Y %*% t(P)
		sigz <- apply(Z, 2, scalefn)
		lambda <- diag(sigz^2)

		list(P=P, lambda=lambda)
	}

	h<-quanf(n=n,p=p,alpha=0.5)
	P<-ogkscatter(x0,scale_est,only.P=TRUE)
    	hsets<-initset(x0,scale_est,P=P,h=h)
	
  	xk<-x0[hsets,,drop=FALSE]
	svd<-classPC(xk,signflip=FALSE) # [P,T,L,r,centerX,meanvct] = classSVD(xk)
        score<-(x0-rep(svd$center,each=n))%*%svd$loadings
        ord<-order(mahalanobisD(score,FALSE,sqrt(abs(svd$eigenvalues))))

  	z<-doScale(x0,center=median,scale=scale_est)
    	z.center<-z$center
    	z.scale<-z$scale
    	z<-z$x
	maxcsteps<-100
 	for(j in 1:maxcsteps){
            score<-(z-rep(svd$center,each=n))%*%svd$loadings
            mah<-mahalanobisD(score,center=FALSE,sd=sqrt(abs(svd$eigenvalues)))
            obs_in_set<-sort.list(mah)[1:h] 
	    xk<-z[obs_in_set,,drop=FALSE]
	    svd<-classPC(xk,signflip=FALSE) # [P,T,L,r,centerX,meanvct] = classSVD(xk)
        }
	return(obs_in_set)
}
