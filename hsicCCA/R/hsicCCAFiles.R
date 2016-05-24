hsicCCAfunc <- function(x,y,Wx=NULL,Wy=NULL,sigmax,sigmay,numiter=20,reltolstop=1e-4) {

  mydistsq <- function(x) {
    xy <- x%*%t(x)
    xx <- diag(xy)
    d <- -2*xy+xx
    d <- t(d)+xx
    d
  }
  numDimX <- ncol(x)
  numDimY <- ncol(y)
  numData <- nrow(x)

  if (is.null(Wx)) Wx <- rnorm(numDimX,sd=0.01)
  if (is.null(Wy)) Wy <- rnorm(numDimY,sd=0.01)
  Wx <- Wx/sqrt(sum(Wx^2))
  Wy <- Wy/sqrt(sum(Wy^2))

  cost <- rep(NA,numiter)
  iter <- 1
  reltol <- Inf

  epsilon <- runif(2,-pi,pi)

    while(iter <= numiter & reltol > reltolstop) {
      projx <- x%*%Wx
      Kx <- exp(-sigmax*mydistsq(projx))
      Kxc <- t(t(Kx)-colMeans(Kx))
      Kxc <- Kxc-rowMeans(Kxc)

      projy <- y%*%Wy
      Ky <- exp(-sigmay*mydistsq(projy))
      Kyc <- t(t(Ky)-colMeans(Ky))
      Kyc <- Kyc-rowMeans(Kyc)

      cost[iter] <-  -sum(Kx*Kyc)/(numData^2)
      cat('iter=',iter,'cost=',cost[iter],'epsilon=',epsilon,'reltol=',reltol,'\n')

      wtSet <- Kx*Kyc
      gradx <- 2*sigmax*t(Wx)%*%sumWtDiff(wtSet,x)/(numData^2)
      wtSet <- Kxc*Ky
      grady <- 2*sigmay*t(Wy)%*%sumWtDiff(wtSet,y)/(numData^2)

      crit <- function(epsilon) {
        hx <- t(gradx)-drop(gradx%*%Wx)*Wx
        nx <- hx/sqrt(sum(hx^2))
        Wx <- Wx*sin(epsilon[1])+nx*cos(epsilon[1])
        hy <- t(grady)-drop(grady%*%Wy)*Wy
        ny <- hy/sqrt(sum(hy^2))
        Wy <- Wy*sin(epsilon[2])+ny*cos(epsilon[2])

        projx <- x%*%Wx
        Kx <- exp(-sigmax*mydistsq(projx))
        Kxc <- t(t(Kx)-colMeans(Kx))
        Kxc <- Kxc-rowMeans(Kxc)
        projy <- y%*%Wy
        Ky <- exp(-sigmay*mydistsq(projy))
        Kyc <- t(t(Ky)-colMeans(Ky))
        Kyc <- Kyc-rowMeans(Kyc)

        -sum(Kx*Kyc)/(numData^2)
      }

      trial <- 0
      cat('running Nelder-Mead. \n')
      epsilonobj <- optim(epsilon,crit,method=c("Nelder-Mead"))
      epsilon <- epsilonobj$par
      reltol <- (epsilonobj$value-cost[iter])/cost[iter]
      while (reltol < 0 & trial < 10) {
      cat('retrying Nelder-Mead. \n')
        trial <- trial+1
        epsilon <- runif(2,-pi,pi)
        epsilonobj <- optim(epsilon,crit,method=c("Nelder-Mead"))
        epsilon <- epsilonobj$par
        reltol <- (epsilonobj$value-cost[iter])/cost[iter]
        }
      if (trial == 10) {
        warning("Nelder-Mead cannot improve the cost. \n")
        break
        }
      iter <- iter+1
      hx <- t(gradx)-drop(gradx%*%Wx)*Wx
      nx <- hx/sqrt(sum(hx^2))
      Wx <- Wx*sin(epsilon[1])+nx*cos(epsilon[1])
      hy <- t(grady)-drop(grady%*%Wy)*Wy
      ny <- hy/sqrt(sum(hy^2))
      Wy <- Wy*sin(epsilon[2])+ny*cos(epsilon[2])
      }
    if (reltol>=0) cost[iter] <- epsilonobj$value
    cat('iter=',iter,'cost=',cost[iter],'epsilon=',epsilon,'reltol=',reltol,'\n')
    return(ls=list(Wx=Wx,Wy=Wy,cost=cost[!is.na(cost)]))
    }



hsicCCA <- function(x,y,M,sigmax=NULL,sigmay=NULL,numrepeat=5,numiter=100,reltolstop=1e-4) {
  Wx <- Wy <- c()
  projx <- diag(1,ncol(x))
  projy <- diag(1,ncol(y))
  for (K in 1:M) {
    fitList <- vector('list',numrepeat)
    if (is.null(sigmax)) sigmax <- median(dist(x))
    if (is.null(sigmay)) sigmay <- median(dist(y))
    for (i in 1:numrepeat) fitList[[i]] <- hsicCCAfunc(x,y,sigmax=sigmax,sigmay=sigmay,reltolstop=reltolstop,numiter=numiter)
    fit <- fitList[[which.min(sapply(fitList,function(x) x$cost[length(x$cost)]))]]
    Wx <- cbind(Wx,projx%*%drop(fit$Wx))
    Wy <- cbind(Wy,projy%*%drop(fit$Wy))
    newprojx <- svd(diag(1,ncol(x))-fit$Wx%*%t(fit$Wx))$u[,1:(ncol(x)-1)]
    x <- x%*%newprojx
    projx <- projx%*%newprojx
    newprojy <- svd(diag(1,ncol(y))-fit$Wy%*%t(fit$Wy))$u[,1:(ncol(y)-1)]
    y <- y%*%newprojy
    projy <- projy%*%newprojy
  }
  return(ls=list(Wx=Wx,Wy=Wy))
}


ktaCCAfunc <- function(x,y,Wx=NULL,Wy=NULL,sigmax,sigmay,numiter=20,reltolstop=1e-4) {

  mydistsq <- function(x) {
    xy <- x%*%t(x)
    xx <- diag(xy)
    d <- -2*xy+xx
    d <- t(d)+xx
    d
  }

    numDimX <- ncol(x)
    numDimY <- ncol(y)
    numData <- nrow(x)

    if (is.null(Wx)) Wx <- rnorm(numDimX,sd=0.01)
    if (is.null(Wy)) Wy <- rnorm(numDimY,sd=0.01)
    Wx <- Wx/sqrt(sum(Wx^2))
    Wy <- Wy/sqrt(sum(Wy^2))

    cost <- rep(NA,numiter)
    iter <- 1
    reltol <- Inf

    epsilon <- runif(2,-pi,pi)

    while(iter <= numiter & reltol > reltolstop) {
      projx <- x%*%Wx
      Kx <- exp(-sigmax*mydistsq(projx))
      Kxc <- t(t(Kx)-colMeans(Kx))
      Kxc <- Kxc-rowMeans(Kxc)

      projy <- y%*%Wy
      Ky <- exp(-sigmay*mydistsq(projy))
      Kyc <- t(t(Ky)-colMeans(Ky))
      Kyc <- Kyc-rowMeans(Kyc)

      numer <- sum(Kx*Kyc)
      denomx <- sum(Kx*Kxc)
      denomy <- sum(Ky*Kyc)
      cost[iter] <-  -log(numer/sqrt(denomx*denomy))
      cat('iter=',iter,'cost=',cost[iter],'epsilon=',epsilon,'reltol=',reltol,'\n')

      wtSet <- Kx*Kyc/numer - Kx*Kxc/denomx
      gradx <- 2*sigmax*t(Wx)%*%sumWtDiff(wtSet,x)
      wtSet <- Kxc*Ky/numer - Ky*Kyc/denomy
      grady <- 2*sigmay*t(Wy)%*%sumWtDiff(wtSet,y)

      crit <- function(epsilon) {
        hx <- t(gradx)-drop(gradx%*%Wx)*Wx
        nx <- hx/sqrt(sum(hx^2))
        Wx <- Wx*sin(epsilon[1])+nx*cos(epsilon[1])
        hy <- t(grady)-drop(grady%*%Wy)*Wy
        ny <- hy/sqrt(sum(hy^2))
        Wy <- Wy*sin(epsilon[2])+ny*cos(epsilon[2])

        projx <- x%*%Wx
        Kx <- exp(-sigmax*mydistsq(projx))
        Kxc <- t(t(Kx)-colMeans(Kx))
        Kxc <- Kxc-rowMeans(Kxc)
        projy <- y%*%Wy
        Ky <- exp(-sigmay*mydistsq(projy))
        Kyc <- t(t(Ky)-colMeans(Ky))
        Kyc <- Kyc-rowMeans(Kyc)

        numer <- sum(Kx*Kyc)
        denomx <- sum(Kx*Kxc)
        denomy <- sum(Ky*Kyc)
        -log(numer/sqrt(denomx*denomy))
      }

      trial <- 0
      cat('running Nelder-Mead. \n')
      epsilonobj <- optim(epsilon,crit,method=c("Nelder-Mead"))
      epsilon <- epsilonobj$par
      reltol <- -(epsilonobj$value-cost[iter])/cost[iter]
      while (reltol < 0 & trial < 10) {
      cat('retrying Nelder-Mead. \n')
        trial <- trial+1
        epsilon <- runif(2,-pi,pi)
        epsilonobj <- optim(epsilon,crit,method=c("Nelder-Mead"))
        epsilon <- epsilonobj$par
        reltol <- -(epsilonobj$value-cost[iter])/cost[iter]
        }
      if (trial == 10) {
        warning("Nelder-Mead cannot improve the cost. \n")
        break
        }
      iter <- iter+1
      hx <- t(gradx)-drop(gradx%*%Wx)*Wx
      nx <- hx/sqrt(sum(hx^2))
      Wx <- Wx*sin(epsilon[1])+nx*cos(epsilon[1])
      hy <- t(grady)-drop(grady%*%Wy)*Wy
      ny <- hy/sqrt(sum(hy^2))
      Wy <- Wy*sin(epsilon[2])+ny*cos(epsilon[2])
      }
    if (reltol>=0) cost[iter] <- epsilonobj$value
    cat('iter=',iter,'cost=',cost[iter],'epsilon=',epsilon,'reltol=',reltol,'\n')
    return(ls=list(Wx=Wx,Wy=Wy,cost=cost[!is.na(cost)]))
    }



ktaCCA <- function(x,y,M,sigmax=NULL,sigmay=NULL,numrepeat=5,numiter=100,reltolstop=1e-4) {
  Wx <- Wy <- c()
  projx <- diag(1,ncol(x))
  projy <- diag(1,ncol(y))
  for (K in 1:M) {
    fitList <- vector('list',numrepeat)
    if (is.null(sigmax)) sigmax <- median(dist(x))
    if (is.null(sigmay)) sigmay <- median(dist(y))
    for (i in 1:numrepeat) fitList[[i]] <- ktaCCAfunc(x,y,sigmax=sigmax,sigmay=sigmay,reltolstop=reltolstop,numiter=numiter)
    fit <- fitList[[which.min(sapply(fitList,function(x) x$cost[length(x$cost)]))]]
    Wx <- cbind(Wx,projx%*%drop(fit$Wx))
    Wy <- cbind(Wy,projy%*%drop(fit$Wy))
    newprojx <- svd(diag(1,ncol(x))-fit$Wx%*%t(fit$Wx))$u[,1:(ncol(x)-1)]
    x <- x%*%newprojx
    projx <- projx%*%newprojx
    newprojy <- svd(diag(1,ncol(y))-fit$Wy%*%t(fit$Wy))$u[,1:(ncol(y)-1)]
    y <- y%*%newprojy
    projy <- projy%*%newprojy
  }
  return(ls=list(Wx=Wx,Wy=Wy))
}

sumWtDiff <- function(Wt,x) {
  numdata <- nrow(x)
  addMat <- Wt
  addMat[lower.tri(addMat,diag=T)] <- 0
  xx <- t(rowSums(addMat+t(addMat))*x)%*%x
  addMat <- addMat[-numdata,]
  crossterm <- t(x[-numdata,])%*%addMat%*%x
  result <- xx-crossterm-t(crossterm)

  addMat <- Wt
  addMat[upper.tri(addMat,diag=T)] <- 0
  addMat <- addMat+t(addMat)
  xx <- t(rowSums(addMat)*x)%*%x
  addMat[lower.tri(addMat,diag=T)] <- 0
  addMat <- addMat[-numdata,]
  crossterm <- t(x[-numdata,])%*%addMat%*%x

  result + xx-crossterm-t(crossterm)
  }
