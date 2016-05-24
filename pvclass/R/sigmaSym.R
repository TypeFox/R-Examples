sigmaSym <- function(X,Y,L,dimension,n,nvec,B=NULL,nu=0,delta=10^(-7),
	prewhitened=FALSE,
	steps=TRUE,nmax=500)
# Computes symmetrized M-estimator of scatter for several classes.
{
	Xnames <- colnames(X)
	Xs <- as.matrix(X)
	if (n <= dimension)
	{
		return("Too few observations!")
	}

        if(is.null(B))
          {
            if (!prewhitened)
              {
		# Prewhiten the data:
		S0 <- 2*cov(Xs)
		B <- SqrtS(S0)
		Xs <- t(qr.solve(B,t(Xs)))
              }
            else
              {
		B <- diag(rep(1,dimension))
              }
        
            # Phase 1: Use a random subset of n
            # observation pairs to get a good
            # starting value:
            if (steps)
              {
		print("Phase 1 (using a random subset of n pairs):",quote=FALSE)
              }
            X0 <- NULL
            for(b in seq_len(L)){
              Xtmp <- unique(X[Y==b, , drop = FALSE])
              Pi <- sample(NROW(Xtmp))
              X0 <- rbind(X0, Xtmp[Pi, , drop = FALSE]
                              - Xtmp[Pi[c(2:NROW(Xtmp),1)], , drop = FALSE])
            }
            
            C <- MVTMLE0(X0,nu,delta,
                         steps=steps,prewhitened=TRUE)$B
            B <- B %*% C
            Xs <- t(qr.solve(C,t(Xs)))
          }
        
	# Phase 2:
	if (steps)
	{
		print("Phase 2 (using all N pairs):",quote=FALSE)
	}
	if (n <= nmax)
	# Use a big data matrix with
	# N observation pairs:
	{
          Xtmp <- lapply(X = seq_len(L),
                         FUN = function(b, mat, Y) unique(mat[which(Y == b), , drop = FALSE]),
                         mat = Xs,
                         Y = Y)
          
          nvec <- sapply(Xtmp, FUN = NROW)
          
          w <- rep(2 / nvec / (n - L), choose(nvec,2))
          XXs <- NULL
          for(b in seq_len(L)) {
            IJ <- combn(nvec[b], 2)
            XXs <- rbind(XXs, Xtmp[[b]][IJ[1,], , drop = FALSE] - Xtmp[[b]][IJ[2,], , drop = FALSE])
          }
          tmp <- MVTMLE0w(XXs,w=w,nu,steps=steps,
                         prewhitened=TRUE)
          B <- B %*% tmp$B
          iter <- tmp$iter
	}
        else
	# Avoid the storage of N observation
	# pairs and perform MVTMLE0 "in loops":
	{
                w <- 2 / nvec / (n - L)
                Psi <- 0
                for(b in seq_len(L)) {
                  Psi <- Psi + LocalPsi(Xs[Y==b, , drop = FALSE],nu,n = nvec[b],dimension,N= 1/w[b])
                }
		Tmp <- eigen(Psi,symmetric=TRUE)
		nG <- sqrt(sum((1 - Tmp$values)^2))
                iter <- 0
		while (nG > delta && iter < 1000)
		{
			iter <- iter + 1
			B <- B %*% Tmp$vectors
			Xs <- Xs %*% Tmp$vectors
                        Ht <- 0
                        for(b in seq_len(L)) {
                          Ht <- Ht + LocalHt(Xs[Y==b,],nu,n = nvec[b],dimension,N= 1/w[b])
                        }

                        a <- qr.solve(diag(Tmp$values) - Ht + (nu == 0)/dimension,Tmp$values - 1)

			# Check whether the new matrix
			# parameter is really better:
			a.half <- a/2
			Xs.new <- t(t(Xs)*exp(-a.half))
                        DL <- 0
                        for(b in seq_len(L)) {
                          DL <- DL + LocalDL(Xs[Y==b,],nu,n = nvec[b],dimension,
                                             N= 1/w[b],Xs.new[Y==b,],a=0)
                        }
                        DL <- DL + sum(a)
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
			} else
			# perform a step of the classical
			# fixed-point algorithm:
			{
				sqrtTmpValues <- sqrt(Tmp$values)
				B <- t(t(B) * sqrtTmpValues)
				Xs <- t(t(Xs) / sqrtTmpValues)
			}

                        Psi <- 0
                        for(b in seq_len(L)) {
                          Psi <- Psi + LocalPsi(Xs[Y==b,],nu,n = nvec[b],dimension,N= 1/w[b])
                        }
			Tmp <- eigen(Psi,symmetric=TRUE) 
			nG <- sqrt(sum((1 - Tmp$values)^2))
                      }
		if (steps)
                  {
			print(cbind(After.iteration=iter,
                                    Norm.of.gradient=nG))
                  }
	}
	S <- tcrossprod(B)
 	if (nu == 0)
 	{
 		S <- S / det(S)^(1/dimension)
 	}
	rownames(S) <- colnames(S) <- Xnames
	return(list(B=B,S=S,iter=iter))
}




LocalPsi <- function(X,nu,n,p,N)
{
	Y <- X[n, , drop = FALSE] - X[n-1, , drop = FALSE]
	Z <- Y/sqrt(nu + sum(Y^2))
	Psi <- tcrossprod(Z)
	for (i in 1:(n-2))
	{
		jj <- (i+1):n
		Y <- t(t(X[jj, , drop = FALSE]) - X[i,])
		denom <- nu + rowSums(Y^2)
		Z <- Y/sqrt(denom)
		Psi <- Psi + crossprod(Z)
	}
	Psi <- Psi * ((nu + p)/N)
	return(Psi)
}


LocalHt <- function(X,nu,n,p,N,evs)
{
	Y <- X[n, , drop = FALSE] - X[n-1, , drop = FALSE]
	Z <- Y^2/(nu + sum(Y^2))
	Ht <- tcrossprod(Z)
	for (i in 1:(n-2))
	{
		jj <- (i+1):n
		Y <- t(t(X[jj, , drop = FALSE]) - X[i,])
		denom <- nu + rowSums(Y^2)
		Z <- Y^2/denom
		Ht <- Ht + crossprod(Z)
	}
        Ht <- Ht * ((nu + p)/N)
	return(Ht)
}


LocalDL <- function(X,nu,n,p,N,Xnew,a)
{
	Y    <- X[n, , drop = FALSE]    - X[n-1, , drop = FALSE]
	Ynew <- Xnew[n,] - Xnew[n-1,]
	DL <- log((nu + sum(Ynew^2))/
	          (nu + sum(Y^2)))
	for (i in 1:(n-2))
	{
		jj <- (i+1):n
		Y    <- t(t(X[jj,])    - X[i,])
		z    <- nu + rowSums(Y^2)
		Ynew <- t(t(Xnew[jj,]) - Xnew[i,])
		znew <- nu + rowSums(Ynew^2)
		DL <- DL + sum(log(znew/z))
	}
	DL <- DL * ((nu + p)/N) + sum(a)
	return(DL)
}

MVTMLE0w <- function(X,w=NULL,nu=0,delta=10^(-7),
	prewhitened=FALSE,
	steps=FALSE)
{
	Xnames <- colnames(X)
	Xs <- as.matrix(X)
	n <- dim(Xs)[1]
	p <- dim(Xs)[2]
	if (length(w) != n)
	{
          w <- rep(1/n,n)
	} else
	{
		w <- w / sum(w)
	}
	if (!prewhitened)
	{
		S0 <- crossprod(Xs)/n
		B <- SqrtS(S0)
		Xs <- t(qr.solve(B,t(Xs)))
	} else
	{
		B <- diag(rep(1,p))
	}

	denom <- nu + rowSums(Xs^2)
	Zs <- Xs / sqrt(denom)
	Psi <- (nu + p) * crossprod(Zs * sqrt(w))
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
			(nu + p) * crossprod(Zs * sqrt(w)) + (nu == 0)/p
		a <- qr.solve(Ht,Tmp$values - 1)

		# Check whether the new matrix
		# parameter is really better:
        a.half <- a/2
		Xs.new <- t(t(Xs)*exp(-a.half))
		denom.new <- nu + rowSums(Xs.new^2)
		DL <- (nu + p)*(log(denom.new/denom) %*% w) +
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
		Zs <- Xs / sqrt(denom)
		Psi <- (nu + p) * crossprod(Zs * sqrt(w))
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

# Auxiliary function:

SqrtS <- function(S)
# For a symmetric, positive definite matrix S,
# this procedure returns a matrix B such that
#    S = B %*% t(B).
{
	res <- eigen(S,symmetric=TRUE)
	B <- t(t(res$vectors) * sqrt(res$values))
	return(B)
}


# Algorithm for the scatter-only problem:

MVTMLE0 <- function(X,nu=0,delta=10^(-7),
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
		B <- SqrtS(S0)
		Xs <- t(qr.solve(B,t(Xs)))
	} else
	{
		B <- diag(rep(1,p))
	}

	denom <- nu + rowSums(Xs^2)
	Zs <- Xs / sqrt(denom)
	Psi <- (nu + p) * crossprod(Zs) / n
	Tmp <- eigen(Psi,symmetric=TRUE)
	nG <- sqrt(sum((1 - Tmp$values)^2))
	iter <- 0
	while (nG > delta && iter < 1000)
	{
		iter <- iter + 1
		B <- B %*% Tmp$vectors
		Xs <- Xs %*% Tmp$vectors
		Zs <- Xs^2/denom
		Ht <- diag(as.matrix(Tmp$values)) -
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
		Zs <- Xs / sqrt(denom)
		Psi <- (nu + p) * crossprod(Zs) / n
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

