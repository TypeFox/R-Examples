stEM <- function(X, V, tol=1e-4, maxiters=500, maxiters.S=1000, print.it=FALSE){
	## check dat
	if( missing(X) ) stop("'X' missing")
	if( missing(V) ) stop("'V' missing")
	
	if(is.data.frame(X) | is.matrix(X))
		X <- data.matrix(X)
	else stop("Data matrix must be of class matrix or data.frame")
	n <- nrow(X)
	p <- ncol(X)

	if(is.data.frame(V) | is.matrix(V))
		V <- data.matrix(V)
	
	wtrim = which(apply(V,1,sum)==0) 
	epsilon = length(wtrim)/n
	
	## init 
	m <- matrix(NA,p,1)
	Sigma <- matrix(0,p,p)
	D <- matrix(NA,n,p)
	Dd <- matrix(NA,n,1)
	Ddtmp <- matrix(NA,n,1)
	autovalues <- matrix(NA,p,1)
	U <- Sigma
  
	## init values for Location and scatter
	Xt <- X					
	Xt[V==0] <- NA 	
	m[,1] <- colMeans(Xt,na.rm=T)
	Sigma <- cov(Xt,use="pairwise.complete.obs")
	s <- eigen(Sigma)
	U <- s$vectors
	autovalues[,1] <- s$values

	Sigma <- U%*%diag(autovalues[,1])%*%t(U)

	## likelihood value
	lik <- sum(ldmvnorm(Xt,m,Sigma),na.rm=TRUE)

	## Start CES step
	likold <- lik - 2*tol
	ii <- 0
	while(lik - likold > tol & ii < maxiters){
		ii <- ii + 1

		iter <- 0
		flag <- FALSE
		while(iter < maxiters.S & flag==FALSE) {
			iter <- iter + 1
			s1 <- sample(which(V==1),1)
			s2 <- sample(which(V[-wtrim,]==0),1)
			Vc <- V
			Vc[s1] <- 0
			Vc[-wtrim,][s2] <- 1

			Xt <- X
			Xt[Vc==0] <- NA
                    	likcand <- sum(ldmvnorm(Xt,m,Sigma),na.rm=TRUE)

			if(likcand > lik){
				V <- Vc
				flag <- TRUE
			}
		}

            # concentration step to update the trim indicator ## 
                       
                	Xt <- X
			Xt[-wtrim,][V[-wtrim,]==0] <- NA
			liks <- ldmvnorm(Xt,m,Sigma)
                        V[wtrim,]=1
                        V[which(liks<quantile(liks,epsilon)),]=0 
                        wtrim=which(liks<quantile(liks,epsilon))

		Dd[,1] <- rep(1,n)
		Dd[which(rowMeans(V) == 0),] <- 0

		## M step
		Xt <- X
		Xt[V==0] <- NA
		XtDd <- sweep(Xt, 1, Dd, "*") 
		VDd <- sweep(V, 1, Dd, "*")
		m[,1] <- colSums(XtDd, na.rm=T)
		m[,1] <- m[,1]/colSums(VDd, na.rm=T)

		Stmp <- matrix(NA, p,p)
		for(h in 1:(p-1)) {
			for(l in (h+1):p) {
	    			Stmp[h,l] <- sum(Dd[,1]*(Xt[,h]-m[h])*(Xt[,l]-m[l]),na.rm=T)
	    			Stmp[h,l] <- Stmp[h,l]/sum(Dd[,1]*V[,h]*V[,l])
				Stmp[l,h] <- Stmp[h,l]
			}
		}
		for(h in 1:p) Stmp[h,h] <- sum(Dd[,1]*(Xt[,h]-m[h])^2,na.rm=T)/sum(Dd[,1]*V[,h])

		## eigenvalue constraint
		s <- eigen(Stmp)
		U <- s$vectors
		autovalues[,1] <- s$values
		#autovalues[,1] <- eigenConst( autovalues[,1], p, restr.fact)
		Sigma <- U%*%diag(autovalues[,1])%*%t(U)
	
		likold <- lik
		lik <- sum(ldmvnorm(Xt, m, Sigma),na.rm=TRUE)

		if( print.it )cat("iter", ii, "; current lik:", lik, "; change in lik:",lik-likold, "\n")
	}
	return(list(mu=m,S=Sigma,V=V,lik=lik,iter=ii))
}


