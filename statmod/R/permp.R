permp <- function(x,nperm,n1,n2,total.nperm=NULL,method="auto",twosided=TRUE)
#	 Exact permutation p-values
#	 Belinda Phipson and Gordon Smyth
#	 16 February 2010. Last modified 27 May 2010.
{
	if(any(x<0)) stop("negative x values")
	if(any(x>nperm)) stop("x cannot exceed nperm")

	if(is.null(total.nperm)) {
		total.nperm <- choose((n1+n2),n1)
		if(n1==n2 & twosided==TRUE) total.nperm <- total.nperm/2
	}

	method <- match.arg(method,c("auto","exact","approximate"))
	if(method=="auto") method <- ifelse(total.nperm>10000,"approximate","exact")

#	exact p-value by summation
	if(method=="exact") {
		p <- (1:total.nperm)/total.nperm
		prob <- rep(p,length(x))
		x2 <- rep(x,each=total.nperm)
		Y <- matrix(pbinom(x2,prob=prob,size=nperm),total.nperm,length(x))
		x[] <- colSums(Y)/total.nperm
	}

#	integral approximation
	else {
		z <- gauss.quad.prob(128,l=0,u=0.5/total.nperm)
		prob <- rep(z$nodes,length(x))
		x2 <- rep(x,each=128)
		Y <- matrix(pbinom(x2,prob=prob,size=nperm),128,length(x))
		int <- 0.5/total.nperm*colSums(z$weights*Y)
		x[] <- (x+1)/(nperm+1)-int
	}
	x
}

