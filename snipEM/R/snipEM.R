# There seems to be a bug in the tclust package. Thus, I switch to a relatively 
# simpler search for the right eigenvalues. See eigenConst. In general, the eigenvalue
# constraints do not seem to play much role on the performance of the estimators. 
# .restr_eigen_C <- function (autovalues, ni.ini, restr.fact=12, zero.tol=1e-16)
# {
	# ##### function parameters:
	# ##### autovalues: matrix containin eigenvalues        
	# ##### ni.ini: vector of current sample sizes for each cluster 
	# ##### restr.fact: eigenvalue restrictions  
	# p <- nrow(autovalues) 
	# K <- ncol(autovalues)                      
	# ret <- .C('RestrictEigenValues', PACKAGE = "tclust"
                # , as.integer (c(p, K))
                # , nParOut = integer (1)
                # , as.double (c (restr.fact, zero.tol))
                # , dParOut = double (1)
                # , EV = as.double (autovalues)
                # , as.double (ni.ini)
      # )
	# if (ret$nParOut[1])
		# matrix (ret$EV, nrow = p)
      # else
		# matrix (0, nrow = p, ncol = K)
# }





snipEM.initialV <- function(X, V, mu0, S0, maxiters.S=100, greedy=TRUE){
	iter <- 0
	flag <- FALSE
	Xt <- X
	Xt[ V == 0 ] <- NA
	lik <- sum(ldmvnorm(Xt, mu0, S0))
	while(iter < maxiters.S & flag==FALSE) {
		iter <- iter + 1
		s1 <- sample(which(V==1),1)
		s2 <- sample(which(V==0),1)
		Vc <- V
		Vc[s1] <- 0
		Vc[s2] <- 1

		Xt <- X
		Xt[Vc==0] <- NA
		likcand <- sum(ldmvnorm(Xt,mu0,S0))

		if(likcand > lik){
			V <- Vc
			## if not a greedy search, stop right after having a V with a larger likelihood
			if(!greedy) flag <- TRUE
		}
	}	
	return(list(V=V, iter=iter))
}



snipEM <- function(X, V, tol=1e-4, maxiters=500, maxiters.S=1000, print.it=FALSE){
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
	if( any(dim(V) != dim(X)) ) stop("'X' and 'V' have non-conforming size")
	epsilon <- sum(V==0)/(n*p)

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
	

	## Impose constraints on the eigenvalues
	#autovalues[,1] <- restr_eigen_C(cbind(autovalues,autovalues), c(n,n), restr.fact)[,1]
	#autovalues[,1] <- eigenConst( autovalues[,1], p, restr.fact)

	Sigma <- U%*%diag(autovalues[,1])%*%t(U)

	## likelihood value
	lik <- sum(ldmvnorm(Xt,m,Sigma))

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
			s2 <- sample(which(V==0),1)
			Vc <- V
			Vc[s1] <- 0
			Vc[s2] <- 1

			Xt <- X
			Xt[Vc==0] <- NA
			likcand <- sum(ldmvnorm(Xt,m,Sigma))

			if(likcand > lik){
				V <- Vc
				flag <- TRUE
			}
		}

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
		lik <- sum(ldmvnorm(Xt, m, Sigma))

		if( print.it )cat("iter", ii, "; current lik:", lik, "; change in lik:",lik-likold, "\n")
	}
	return(list(mu=m,S=Sigma,V=V,lik=lik,iter=ii))
}
