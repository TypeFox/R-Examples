#  NUMERICAL INTEGRATION

gauss.quad <- function(n,kind="legendre",alpha=0,beta=0)
#	Calculate nodes and weights for Gaussian quadrature.
#	Adapted from Netlib routine gaussq.f
#	Gordon Smyth, Walter and Eliza Hall Institute
#	Suggestion from Stephane Laurent 6 Aug 2012
#	4 Sept 2002. Last modified 7 Aug 2012.
{
	n <- as.integer(n)
	if(n<0) stop("need non-negative number of nodes")
	if(n==0) return(list(nodes=numeric(0), weights=numeric(0)))
	kind <- match.arg(kind,c("legendre","chebyshev1","chebyshev2","hermite","jacobi","laguerre"))
	i <- 1:n
	i1 <- i[-n] # 1:(n-1)
	switch(kind, legendre={
		lnmuzero <- log(2)
		a <- rep(0,n)
		b <- i1/sqrt(4*i1^2-1)
	}, chebyshev1={
		lnmuzero <- log(pi)
		a <- rep(0,n)
		b <- rep(0.5,n-1)
		b[1] <- sqrt(0.5)
	}, chebyshev2={
		lnmuzero <- log(pi/2)
		a <- rep(0,n)
		b <- rep(0.5,n-1)
	}, hermite={
		lnmuzero <- log(pi)/2
		a <- rep(0,n)
		b <- sqrt(i1/2)
	}, jacobi={
		ab <- alpha+beta
#		muzero <- 2^(ab+1) * gamma(alpha+1) * gamma(beta+1) / gamma(ab+2)
		lnmuzero <- (ab+1)*log(2) + lgamma(alpha+1) + lgamma(beta+1) - lgamma(ab+2)
		a <- i
		a[1] <- (beta-alpha)/(ab+2)
		i2 <- 2:n
		abi <- ab+2*i2
		a[i2] <- (beta^2-alpha^2)/(abi-2)/abi
		b <- i1
		b[1] <- sqrt(4*(alpha+1)*(beta+1)/(ab+2)^2/(ab+3))
		i2 <- i1[-1] # 2:(n-1)
		abi <- ab+2*i2
		b[i2] <- sqrt(4*i2*(i2+alpha)*(i2+beta)*(i2+ab)/(abi^2-1)/abi^2)
	}, laguerre={
		a <- 2*i-1+alpha
		b <- sqrt(i1*(i1+alpha))
		lnmuzero <- lgamma(alpha+1)
	})
	b <- c(b,0)
	z <- rep.int(0,n)
	z[1] <- 1
	ierr <- 0L
   out <- .Fortran("gausq2",n,as.double(a),as.double(b),as.double(z),ierr,PACKAGE="statmod")
	x <- out[[2]]
	w <- out[[4]]
	w <- exp(lnmuzero + 2*log(abs(w)))
	list(nodes=x,weights=w)
}

gauss.quad.prob <- function(n,dist="uniform",l=0,u=1,mu=0,sigma=1,alpha=1,beta=1)
#	Calculate nodes and weights for Guassian quadrature using probability densities.
#	Adapted from Netlib routine gaussq.f
#	Gordon Smyth, Walter and Eliza Hall Institute
#	Corrections for n=1 and n=2 by Spencer Graves, 28 Dec 2005
#	4 Sept 2002. Last modified 7 Aug 2012.
{
	n <- as.integer(n)
	if(n<0) stop("need non-negative number of nodes")
	if(n==0) return(list(nodes=numeric(0), weights=numeric(0)))
	dist <- match.arg(dist,c("uniform","beta1","beta2","normal","beta","gamma"))
	if(n==1){
		switch(dist,
			uniform={x <- (l+u)/2},
			beta1=,beta2=,beta={x <- alpha/(alpha+beta)},
			normal={x <- mu},
			gamma={x <- alpha*beta}
		)
		return(list(nodes=x, weights=1))
	}
	if(dist=="beta" && alpha==0.5 && beta==0.5) dist <- "beta1"
	if(dist=="beta" && alpha==1.5 && beta==1.5) dist <- "beta2"
	i <- 1:n
	i1 <- 1:(n-1)
	switch(dist, uniform={
		a <- rep(0,n)
		b <- i1/sqrt(4*i1^2-1)
	}, beta1={
		a <- rep(0,n)
		b <- rep(0.5,n-1)
		b[1] <- sqrt(0.5)
	}, beta2={
		a <- rep(0,n)
		b <- rep(0.5,n-1)
	}, normal={
		a <- rep(0,n)
		b <- sqrt(i1/2)
	}, beta={
		ab <- alpha+beta
		a <- i
		a[1] <- (alpha-beta)/ab
		i2 <- 2:n
		abi <- ab-2+2*i2
		a[i2] <- ((alpha-1)^2-(beta-1)^2)/(abi-2)/abi
		b <- i1
		b[1] <- sqrt(4*alpha*beta/ab^2/(ab+1))
		i2 <- i1[-1] # 2:(n-1)
		abi <- ab-2+2*i2
		b[i2] <- sqrt(4*i2*(i2+alpha-1)*(i2+beta-1)*(i2+ab-2)/(abi^2-1)/abi^2)
	}, gamma={
		a <- 2*i+alpha-2
		b <- sqrt(i1*(i1+alpha-1))
	})
	b <- c(b,0)
	z <- rep.int(0,n)
	z[1] <- 1
	ierr <- 0L
	out <- .Fortran("gausq2",n,as.double(a),as.double(b),as.double(z),ierr,PACKAGE="statmod")
	x <- out[[2]]
	w <- out[[4]]^2
	switch(dist,
		uniform = x <- l+(u-l)*(x+1)/2,
		beta1=,beta2=,beta = x <- (x+1)/2,
		normal = x <- mu + sqrt(2)*sigma*x,
		gamma = x <- beta*x)
	list(nodes=x,weights=w)
}
