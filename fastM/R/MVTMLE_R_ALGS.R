
MVTMLE0r <- function(X,nu=0,delta=1e-06,
	prewhitened=FALSE,
	steps=FALSE)
{
	Xnames <- colnames(X)
	Xs <- as.matrix(X)
	n <- dim(Xs)[1]
	p <- dim(Xs)[2]
	if (!prewhitened)
	{
		S0 <- crossprod(Xs)/n
		tmp <- eigen(S0,symmetric=TRUE)
		tmpv <- sqrt(tmp$values)
		B <- t(t(tmp$vectors) * tmpv)
		Xs <- Xs %*% tmp$vectors
		Xs <- t(t(Xs) / tmpv)
	} else
	{
		B <- diag(rep(1,p))
	}

	denom <- nu + rowSums(Xs^2)
	Ys <- Xs / sqrt(denom)
	Psi <- (nu + p) * crossprod(Ys) / n
	Tmp <- eigen(Psi,symmetric=TRUE)
	nG <- sqrt(sum((1 - Tmp$values)^2))
	iter <- 0
	while (nG > delta && iter < 1000)
	{
		iter <- iter + 1
		B <- B %*% Tmp$vectors
		Xs <- Xs %*% Tmp$vectors
		Zs <- Xs^2/denom
		Ht <- diag(Tmp$values) -
			(nu + p) * crossprod(Zs)/n + (nu == 0)/p
		a <- qr.solve(Ht,Tmp$values - 1)

		# Check whether the new matrix
		# parameter is really better:
        a.half <- a/2
		Xs.new <- t(t(Xs)*exp(-a.half))
		denom.new <- nu + rowSums(Xs.new^2)
		DL <- (nu + p)*mean(log(denom.new/denom)) +
			sum(a)
		DL0 <- sum(a*(1 - Tmp$values))/4
		if (steps)
		{
			print(cbind(Before.iteration=iter,
				Norm.of.gradient=nG,
				Diff.M2LL.with.PN=DL))
		}
		if (DL <= DL0)
		{
			B <- t(t(B) * exp(a.half))
			Xs <- Xs.new
			denom <- denom.new
		} else
		# perform a step of the classical
		# fixed-point algorithm:
		{
			sqrtTmpValues <- sqrt(Tmp$values)
			B <- t(t(B) * sqrtTmpValues)
			Xs <- t(t(Xs) / sqrtTmpValues)
			denom <- nu + rowSums(Xs^2)
		}
		Ys <- Xs / sqrt(denom)
		Psi <- (nu + p) * crossprod(Ys) / n
		Tmp <- eigen(Psi,symmetric=TRUE)
		nG <- sqrt(sum((1 - Tmp$values)^2))
	}
	if (steps)
	{
		print(cbind(After.iteration=iter,
			Norm.of.gradient=nG))
	}
	S <- tcrossprod(B)
	rownames(S) <- colnames(S) <- Xnames
	if (nu == 0)
	{
		S <- S / det(S)^(1/p)
	}
	return(list(B=B,S=S,iter=iter))
}



# Alternative algorithms for the scatter-only
# problem:

MVTMLE0r_FP0 <- function(X,nu=0,delta=1e-06,
	steps=FALSE)
{
    X <- as.matrix(X)
	n <- dim(X)[1]
	p <- dim(X)[2]
	Ip <- diag(rep(1,p))
	S <- crossprod(X)/n

	denom <- nu + rowSums(X*t(qr.solve(S,t(X))))
	Z <- X / sqrt(denom)
	S.new <- (nu + p) * crossprod(Z)/n
	error <- norm(qr.solve(S,S.new)-Ip,type="F")
	iter <- 0
	while (error > delta && iter < 1000)
	{
		iter <- iter+1
		if (steps)
		{
			print(cbind(Before.iteration=iter,
				error=error))
		}
		S <- S.new
		denom <- nu + rowSums(X*t(qr.solve(S,t(X))))
		Z <- X / sqrt(denom)
		S.new <- (nu + p) * crossprod(Z)/n
		error <- norm(qr.solve(S,S.new)-Ip,type="F")
	}
	if (steps)
	{
		print(cbind(After.iteration=iter,error=error))
	}
	S <- S.new
	if (nu == 0)
	{
		S <- S / det(S)^(1/p)
	}
	return(list(S=S,iter=iter))
}

MVTMLE0r_FP <- function(X,nu=0,delta=1e-06,
	steps=FALSE)
{
    X <- as.matrix(X)
	n <- dim(X)[1]
	p <- dim(X)[2]
	S <- crossprod(X)/n
	tmp <- eigen(S,symmetric=TRUE)
	tmpv <- sqrt(tmp$values)
	B <- t(t(tmp$vectors) * tmpv)
	Xs <- X %*% tmp$vectors
	Xs <- t(t(Xs) / tmpv)

	denom <- nu + rowSums(Xs^2)
	Ys <- Xs / sqrt(denom)
	Psi <- (nu + p) * crossprod(Ys)/n
	tmp <- eigen(Psi,symmetric=TRUE)
	nG <- sqrt(sum((1 - tmp$values)^2))
	iter <- 0
	while (nG > delta && iter < 1000)
	{
		iter <- iter+1
		if (steps)
		{
			print(cbind(Before.iteration=iter,
				Norm.of.gradient=nG))
		}

		tmpv <- sqrt(tmp$values)
		B <- B %*% t(t(tmp$vectors)*tmpv)
		Xs <- Xs %*% tmp$vectors
		Xs <- t(t(Xs)/tmpv)

		denom <- nu + rowSums(Xs^2)
		Ys <- Xs/sqrt(denom)
		Psi <- (nu + p) * crossprod(Ys)/n
		tmp <- eigen(Psi,symmetric=TRUE)
		nG <- sqrt(sum((1 - tmp$values)^2))
	}
	if (steps)
	{
		print(cbind(After.iteration=iter,
			Norm.of.gradient=nG))
	}
	S <- tcrossprod(B)
	if (nu == 0)
	{
		S <- S / det(S)^(1/p)
	}
	return(list(S=S,iter=iter))
}

MVTMLE0r_G <- function(X,nu=0,delta=1e-06,
	steps=FALSE)
{
    X <- as.matrix(X)
	n <- dim(X)[1]
	p <- dim(X)[2]
	S <- crossprod(X)/n
	tmp <- eigen(S)
	tmpv <- sqrt(tmp$values)
	B <- t(t(tmp$vectors)*tmpv)
	Xs <- X %*% tmp$vectors
	Xs <- t(t(Xs)/tmpv)
	
	denom <- nu + rowSums(Xs^2)
	Ys <- Xs/sqrt(denom)
	Psi <- (nu + p) * crossprod(Ys)/n
	Tmp <- eigen(Psi)
	Gt <- 1 - Tmp$values
	nGt <- sqrt(sum(Gt^2))
	iter <- 0
	while (nGt > delta && iter < 1000)
	{
		iter <- iter + 1
		Xs <- Xs %*% Tmp$vectors
		B <- B %*% Tmp$vectors
		Zs <- Xs^2/denom

		HGG <- sum(Tmp$values * Gt^2) -
			(nu + p) * sum(colSums(t(Zs)*Gt)^2)/n
		At <- -(nGt^2/HGG)*Gt

		# Check whether the new matrix
		# parameter is really better:
		At.half <- At/2
		Xs.new <- t(t(Xs)*exp(-At.half))
		denom.new <- nu + rowSums(Xs.new^2)
		DL <- (nu + p)*mean(log(denom.new/denom)) +
			sum(At)
		DL0 <- - nGt^2/4
		if (steps)
		{
			print(cbind(Before.iteration=iter,
				Norm.of.gradient=nGt,
				Diff.M2LL.with.G=DL))
		}
		if (DL <= DL0)
		{
			B <- t(t(B)*exp(At.half))
			Xs <- Xs.new
			denom <- denom.new
		} else
		# perform a step of the classical
		# fixed-point algorithm:
		{
			sqrtTmpValues <- sqrt(Tmp$values)
			B <- t(t(B) * sqrtTmpValues)
			Xs <- t(t(Xs) / sqrtTmpValues)
			denom <- nu + rowSums(Xs^2)
		}
		
		Ys <- Xs/sqrt(denom)
		Psi <- (nu + p) * crossprod(Ys)/n
		Tmp <- eigen(Psi)
		Gt <- 1 - Tmp$values
		nGt <- sqrt(sum(Gt^2))
	}
	if (steps)
	{
		print(cbind(After.iteration=iter,
			Norm.of.gradient=nGt))
	}
	S <- tcrossprod(B)
	if (nu == 0)
	{
		S <- S / det(S)^(1/p)
	}
	return(list(S=S,iter=iter))
}

MVTMLE0r_CG <- function(X,nu=0,delta=1e-06,
	steps=FALSE)
{
    X <- as.matrix(X)
	n <- dim(X)[1]
	p <- dim(X)[2]
	S <- crossprod(X)/n
	tmp <- eigen(S)
	tmpv <- sqrt(tmp$values)
	B <- t(t(tmp$vectors)*tmpv)
	Xs <- X %*% tmp$vectors
	Xs <- t(t(Xs)/tmpv)
	
	denom <- nu + rowSums(Xs^2)
	Ys <- Xs/sqrt(denom)
	Psi <- (nu + p) * crossprod(Ys)/n
	A1 <- matrix(0,p,p)
	A2 <- diag(rep(1,p)) - Psi
	error <- norm(A2,type="F")
	iter <- 0
	while (error > delta && iter < 1000)
	{
		if (steps)
		{
			print(cbind(Before.iteration=iter,error=error))
		}
		
		H11 <- (nu + p) *
			sum((rowSums((Ys %*% A1)^2) -
			     rowSums(Ys * (Ys %*% A1))^2))/n
		H12 <- (nu + p) *
			sum((rowSums((Ys %*% A1) * (Ys %*% A2)) -
			     rowSums(Ys * (Ys %*% A1))*
			     rowSums(Ys * (Ys %*% A2))))/n
		H22 <- (nu + p) *
			sum((rowSums((Ys %*% A2)^2) -
			     rowSums(Ys * (Ys %*% A2))^2))/n

		if (H12^2 >= 0.99 * H11*H22)
		{
			A <- (sum(A2^2) / H22) * A2
		}
		else
		{
			dd <- H11*H22 - H12^2
			y1 <- sum(A1*A2)
			y2 <- sum(A2^2)
			t1 <- - (H22*y1 - H12*y2)/dd
			t2 <- - (H11*y2 - H12*y1)/dd
			A <- t1*A1 + t2*A2
		}
		tmp <- eigen(A)
		tmpv <- exp(tmp$values/2)
		B <- B %*% t(t(tmp$vectors)*tmpv)
		Xs <- Xs %*% tmp$vectors
		Xs <- t(t(Xs)/tmpv)

		denom <- nu + rowSums(Xs^2)
		Ys <- Xs/sqrt(denom)
		Psi <- (nu + p) * crossprod(Ys)/n
		A1 <- A
		A2 <- diag(rep(1,p)) - Psi
		error <- norm(A2,type="F")

		iter <- iter + 1
	}
	if (steps)
	{
		print(cbind(iteration=iter,error=error))
	}
	S <- tcrossprod(B)
	if (nu == 0)
	{
		S <- S / det(S)^(1/p)
	}
	return(list(S=S,iter=iter))
}
