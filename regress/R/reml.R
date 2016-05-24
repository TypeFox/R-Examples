## when two matrices are passed to regress this is also called
## to evaluate the REML at certain values of gamma and find a
## good place to start the regress algorithm

reml <- function(lambda, y, X, V0, V1,verbose=0){

  if(is.null(dim(y)))
    {
      isNA <- is.na(y)
      y <- y[isNA==F]
    } else {
      isNA <- apply(y,1,is.na)

      if(is.matrix(isNA))  {
        isNA <- as.logical(apply(isNA,2,sum))
      }
      y <- y[isNA==F,]
    }
  V0 <- V0[isNA==F,isNA==F]
  V1 <- V1[isNA==F,isNA==F]
  X <- X[isNA==F,]
  X <- as.matrix(X)

  qr <- qr(X)
  ##print(qr$rank)
  n <- dim(as.matrix(y))[1]
  In <- diag(1,n)

  X <- matrix(X[, qr$pivot[1:qr$rank]],n,qr$rank)
  llik <- rep(0, length(lambda))
  if(is.null(dim(y))) q <- 1 else q <- dim(y)[2]

  n <- dim(X)[1]
  if(missing(V0)) V0 <- diag(rep(1, n), n, n)
  rank <- n - qr$rank
  ##if(verbose==1) cat("n-p =",n,"-",qr$rank,"=",rank,"\n")
  for(i in 1:length(lambda))
    {
      if(verbose>=2) cat(lambda[i],"\n")

      Sigma <- (1-lambda[i])*V0 + lambda[i] * V1
      ##cholesky <- try(chol(Sigma, pivot=FALSE), silent=TRUE)
      ##if(class(cholesky) == "try-error" || min(diag(cholesky)) < 1e-9) return(1e+32)
      ##W <- chol2inv(cholesky)
      ##WX <- W %*% X
      WX <- solve(Sigma,cbind(X,In))
      W <- WX[,(dim(X)[2]+1):dim(WX)[2]]
      WX <- WX[,1:dim(X)[2]]
      XtWX <- t(X)%*%WX
      WQ <- W - WX%*%solve(XtWX,t(WX))
      rss <- t(y) %*% WQ %*% y
      logdetrss <- sum(log(eigen(rss)$values[1:q]))
      eVals <- eigen(WQ,symmetric=TRUE,only.values=TRUE)$values[1:rank]
      ldet <- sum(log(eVals))
      llik[i] <- Re(ldet*q/2 - rank*logdetrss/2)
    }
  imax <- sort.list(-llik)[1]
  lambdamax <- lambda[imax]
  curv <- 0
  if(imax > 1 && imax < length(lambda)){
    delta <- (lambda[imax+1] - lambda[imax-1])/2
    slope <-  (llik[imax+1] - llik[imax-1])/2
    curv <- llik[imax-1] -2*llik[imax] + llik[imax+1]
    lambdamax <- lambdamax - slope/curv * delta
    curv <- -curv/delta^2
  }
  lamMax <- lambdamax
  Sigma <- (1-lamMax)*V0 + lamMax * V1
  ##cholesky <- try(chol(Sigma, pivot=FALSE), silent=TRUE)
  ##if(class(cholesky) == "try-error" || min(diag(cholesky)) < 1e-9) return(1e+32)
  ##W <- chol2inv(cholesky)
  ##WX <- W %*% X
  WX <- solve(Sigma,cbind(X,In))
  W <- WX[,(dim(X)[2]+1):dim(WX)[2]]
  WX <- WX[,1:dim(X)[2]]
  XtWX <- t(X)%*%WX
  FItWX <- solve(XtWX,t(WX))
  WQ <- W - WX%*%FItWX
  rss <- t(y) %*% WQ %*% y
  beta <- FItWX %*% y

  list(llik=as.numeric(llik),rms=rss/rank, beta=beta, gamma=lambda, gamMax=lambdamax,W=W)
}

remlOptimize <- function(y, X, V0, V1,verbose=0,...){

  if(is.null(dim(y)))
    {
      isNA <- is.na(y)
      y <- y[isNA==F]
    } else {
      isNA <- apply(y,1,is.na)

      if(is.matrix(isNA))  {
        isNA <- as.logical(apply(isNA,2,sum))
      }
      y <- y[isNA==F,]
    }
  V0 <- V0[isNA==F,isNA==F]
  V1 <- V1[isNA==F,isNA==F]
  X <- X[isNA==F,]
  X <- as.matrix(X)

  qr <- qr(X)
  ##print(qr$rank)
  n <- dim(as.matrix(y))[1]
  In <- diag(1,n)

  X <- matrix(X[, qr$pivot[1:qr$rank]],n,qr$rank)
  if(is.null(dim(y))) q <- 1 else q <- dim(y)[2]

  n <- dim(X)[1]
  if(missing(V0)) V0 <- diag(rep(1, n), n, n)
  rank <- n - qr$rank
  ##if(verbose==1) cat("n-p =",n,"-",qr$rank,"=",rank,"\n")

  f <- function(lambda,verbose=verbose) {
    if(verbose>=2) cat(lambda,"\n")
    Sigma <- (1-lambda)*V0 + lambda * V1
    ##cholesky <- try(chol(Sigma, pivot=FALSE), silent=TRUE)
    ##if(class(cholesky) == "try-error" || min(diag(cholesky)) < 1e-9) return(1e+32)
    ##W <- chol2inv(cholesky)
    ##WX <- W %*% X
    WX <- solve(Sigma,cbind(X,In))
    W <- WX[,(dim(X)[2]+1):dim(WX)[2]]
    WX <- WX[,1:dim(X)[2]]
    XtWX <- t(X)%*%WX
    WQ <- W - WX%*%solve(XtWX,t(WX))
    rss <- t(y) %*% WQ %*% y
    logdetrss <- sum(log(eigen(rss)$values[1:q]))
    eVals <- eigen(WQ,symmetric=TRUE,only.values=TRUE)$values[1:rank]
    ldet <- sum(log(eVals))
    llik <- Re(ldet*q/2 - rank*logdetrss/2)
    llik
  }

  res <- optimize(f,interval=c(0,1),maximum=TRUE,verbose=verbose,...)
  lamMax <- res$maximum
  llikMax <- res$objective

  Sigma <- (1-lamMax)*V0 + lamMax * V1
  ##cholesky <- try(chol(Sigma, pivot=FALSE), silent=TRUE)
  ##if(class(cholesky) == "try-error" || min(diag(cholesky)) < 1e-9) return(1e+32)
  ##W <- chol2inv(cholesky)
  ##WX <- W %*% X
  WX <- solve(Sigma,cbind(X,In))
  W <- WX[,(dim(X)[2]+1):dim(WX)[2]]
  WX <- WX[,1:dim(X)[2]]
  XtWX <- t(X)%*%WX
  FItWX <- solve(XtWX,t(WX))
  WQ <- W - WX%*%FItWX
  rss <- t(y) %*% WQ %*% y
  beta <- FItWX %*% y

  list(llik=as.numeric(llikMax),rms=rss/rank, beta=beta, gamMax=lamMax,W=W)
}

