GroupwiseWMeans <- function(X,Y,L,n,W=rep(1,length(Y)))
{
	Id <- diag(rep(1,L))
	D <- Id[Y,]

	sw <- sqrt(W)
	M <- qr.solve(sw*D,sw*X)
	return(M)
}



MVTMLE.LDA <- function(X,Y,L,n,nu=1,
	M0=NULL,B0=NULL, # Potential starting values;
	                 # need both or none.
	delta=10^(-7),steps=FALSE)
{
	Id <- diag(rep(1,dim(X)[2]))
	if (is.null(M0))
	{	# Starting point: B = Id:
		B <- Id
		Xs <- X
		M <- GroupwiseWMeans(Xs,Y,L,n)
	}
	else
	{
		B <- B0
		Xs <- t(qr.solve(B0,t(X)))
		M <- t(qr.solve(B0,t(M0)))
	}
	Xsc <- Xs - M[Y,]
	W <- 1/(nu + rowSums(Xsc^2))
	deltaM <- sqrt(sum(GroupwiseWMeans(Xsc,Y,L,n,W)^2))
	BF <- MVTMLE0(X=Xsc,nu=nu,
		delta=delta,prewhitened=!is.null(M0))$B
        
	deltaS <- sqrt(sum((tcrossprod(BF) - Id)^2))
	iter <- 0
	while (iter < 500 & deltaM + deltaS > delta)
	{
		if (steps)
		{
			print(c(iter=iter,deltaM=deltaM,deltaS=deltaS))
		}
		# Update B:
		B <- B %*% BF
		# Update Xs:
		Xs <- t(qr.solve(BF,t(Xs)))
		# Update M and W:
		M <- GroupwiseWMeans(Xs,Y,L,n,W)
		Xsc <- Xs - M[Y,]
		W <- 1/(nu + rowSums(Xsc^2))
		# Stopping criteria and new
		# potential increment for B:
		deltaM <- sqrt(sum(GroupwiseWMeans(Xsc,Y,L,n,W)^2))
		BF <- MVTMLE0(X=Xsc,nu=nu,
			delta=delta,prewhitened=TRUE)$B
		deltaS <- sqrt(sum((tcrossprod(BF) - Id)^2))
		iter <- iter + 1
	}
	return(list(M=t(B%*%t(M)),B=B,S=tcrossprod(B)))
}
