MVTMLE_symm <- function(X,nu=0,delta=10^(-7),
                        prewhitened=FALSE,
                        steps=TRUE,nmax=500)
  # Computes symmetrized M-estimator of scatter.
{
  Xnames <- colnames(X)
  Xs <- as.matrix(X)
  n <- dim(Xs)[1]
  p <- dim(Xs)[2]
  if (n <= p)
  {
    return("Too few observations!")
  }
  if (!prewhitened)
  {
    # Prewhiten the data:
    S0 <- 2*cov(Xs)
    B <- SqrtS(S0)
    Xs <- t(qr.solve(B,t(Xs)))
  }
  else
  {
    B <- diag(rep(1,p))
    Xs <- X
  }
  
  # Phase 1: Use a random subset of n
  # observation pairs to get a good
  # starting value:
  if (steps)
  {
    print("Phase 1 (using a random subset of n pairs):",quote=FALSE)
  }
  #Pi <- sample(n)
  X0 <- Xs[,] - Xs[c(2:n,1),]
  C <- MVTMLE0r(X0,nu,delta,
                steps=steps,prewhitened=TRUE)$B
  B <- B %*% C
  Xs <- t(qr.solve(C,t(Xs)))
  
  # Phase 2:
  if (steps)
  {
    print("Phase 2 (using all N pairs):",quote=FALSE)
  }
  if (n <= nmax)
    # Use a big data matrix with
    # N observation pairs:
  {
    IJ <- combn(n,2)
    XXs <- Xs[IJ[1,],] - Xs[IJ[2,],]
    tmp <- MVTMLE0r(XXs,nu,delta=delta,
                    prewhitened=TRUE,steps=steps)
    B <- B %*% tmp$B
    iter <- tmp$iter
  }
  else
    # Avoid the storage of N observation
    # pairs and perform MVTMLE0 "in loops":
  {
    N <- n*(n-1)/2
    
    Psi <- LocalPsi(Xs,nu,n,p,N)
    Tmp <- eigen(Psi,symmetric=TRUE)
    nG <- sqrt(sum((1 - Tmp$values)^2))
    iter <- 0
    while (nG > delta && iter < 1000)
    {
      iter <- iter + 1
      B <- B %*% Tmp$vectors
      Xs <- Xs %*% Tmp$vectors
      a <- LocalA(Xs,nu,n,p,N,Tmp$values)
      
      # Check whether the new matrix
      # parameter is really better:
      a.half <- a/2
      Xs.new <- t(t(Xs)*exp(-a.half))
      DL <- LocalDL(Xs,nu,n,p,N,Xs.new,a)
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
      Psi <- LocalPsi(Xs,nu,n,p,N)
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
    S <- S / det(S)^(1/p)
  }
  rownames(S) <- colnames(S) <- Xnames
  return(list(B=B,S=S,iter=iter))
}


LocalPsi <- function(X,nu,n,p,N)
{
	Y <- X[n,] - X[n-1,]
	Z <- Y/sqrt(nu + sum(Y^2))
	Psi <- tcrossprod(Z)
	for (i in 1:(n-2))
	{
		jj <- (i+1):n
		Y <- t(t(X[jj,]) - X[i,])
		denom <- nu + rowSums(Y^2)
		Z <- Y/sqrt(denom)
		Psi <- Psi + crossprod(Z)
	}
	Psi <- Psi * ((nu + p)/N)
	return(Psi)
}

LocalA <- function(X,nu,n,p,N,evs)
{
	Y <- X[n,] - X[n-1,]
	Z <- Y^2/(nu + sum(Y^2))
	Ht <- tcrossprod(Z)
	for (i in 1:(n-2))
	{
		jj <- (i+1):n
		Y <- t(t(X[jj,]) - X[i,])
		denom <- nu + rowSums(Y^2)
		Z <- Y^2/denom
		Ht <- Ht + crossprod(Z)
	}
	Ht <- diag(evs) -
		Ht * ((nu + p)/N) + (nu == 0)/p
	a <- qr.solve(Ht,evs - 1)
	return(a)
}

LocalDL <- function(X,nu,n,p,N,Xnew,a)
{
	Y    <- X[n,]    - X[n-1,]
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