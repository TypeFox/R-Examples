AEnet.aft <-
function(X, Y, delta, weight, lambda2, maxit=10)
{
	n <- nrow(X) # number of samples
	p <- ncol(X) # number of predictors
	if(n != length(delta) || n != length(Y))
	stop("dimensions of X, Y and censorship don't match!")
	weight[weight==0]<-0.001
	w <-1/weight
	kw <- aft.kmweight(Y,delta)$kmwts
	XW <- apply(as.matrix(X[delta == 1,] * kw[delta == 1]), 2,
	sum) / sum(kw[delta == 1])
	YW <- sum(Y[delta == 1] * kw[delta == 1]) /
	sum(kw[delta == 1])
	for(i in 1:n)
	X[i,] <- X[i,] - XW
	X <- as.matrix(sqrt(kw) * X)
	Y <- sqrt(kw) * (Y - YW)
	Xextra<-diag(sqrt(lambda2), p)
	Yextra<-c(rep(0,p))
	delta.extra<-c(rep(1,p))

	X<-rbind(X, Xextra)
	for (j in 1:p)
	X[,j] <-X[,j]/w[j]
	Y<-c(Y, Yextra)
	delta<-c(delta, delta.extra)
	meanx <- apply(X,2,mean)
	normx <- sqrt(apply(X^2,2,sum))
	meany <- mean(Y)
	Y <- Y-meany
	X <- t(t(X)/normx)
    XX <- t(X)%*%X
    h <- rep(0,p)
    beta <- rbind(rep(0,p))
    cor <- t(X)%*%Y
    it <- 1
    potential <- active <- old.active <- rep(F,p)
    NZ <- rep(F,2*p)

get.h <-
		function(X,Y,beta,cor,active,it,XX,NZ){
    p <- ncol(X)
    K <- sum(active)
    XX <- matrix((XX*sign(cor))[active,],nrow=K)
    A.mat <-  cbind(-XX,XX,diag(K))
    NZ <- c(NZ,rep(F,K))
    exclude <- c(NZ[((p+1):(2*p))],NZ[1:p],rep(F,K))
    cor.sign <- sign(cor[which.min(-abs(cor))])
    if (K>1){
      B1.inv <- solve(A.mat[-K,NZ])
      B2 <- A.mat[K,NZ]
      B2.B1inv <- B2%*%B1.inv
      alpha <- as.vector(B2.B1inv%*%rep(1,K-1)-1)
      B1.BK <- B1.inv%*%A.mat[-K,!NZ & !exclude]
      q <- A.mat[K,!NZ& ! exclude]-B2%*%B1.BK
      cost <- c(rep(1,2*p),rep(0,K))
      eval <- alpha*(t(cost[NZ])%*%B1.BK-cost[!NZ & !exclude])/q
      eval[alpha*q<0] <- -10^100
      best <- which.min(-eval)
      new.col <- (1:(2*p+K))[!NZ & !exclude][best]
      old.cols <- (1:(2*p+K))[NZ]
      NZ[new.col] <- T
      q <- q[best]
      B.inv <- B1.inv*q+B1.BK[,best]%*%B2.B1inv
      B.inv <- rbind(cbind(B.inv,-B1.BK[,best]),cbind(-B2.B1inv,1))/q
      B.inv <- B.inv[order(c(old.cols,new.col)),]
    }
    else{
      NZ[which.min(-abs(cor))-p*(cor.sign-1)/2] <- T
      B.inv <- 1/A.mat[,NZ]}
    active.beta <- NZ[1:p] | NZ[(p+1):(2*p)]
    h <- rep(0,p)
    beta.pos <- (1:p)[NZ[1:p]]
    beta.neg <- (1:p)[NZ[(p+1):(2*p)]]
    h[c(beta.pos,beta.neg)] <- -(B.inv%*%rep(1,K))[1:sum(active.beta)]
    h[beta.neg] <- -h[beta.neg]
    slack <- 0
    if (sum(NZ[-(1:(2*p))])>0)
      slack <- (1:p)[active][NZ[-(1:(2*p))]]
    list(h=h,slack=slack,NZ=NZ[1:(2*p)])}

    active[which.max(abs(t(X)%*%Y))] <- T
    while (max(abs(cor))>10^-5 & it<=maxit){
      res <- Y-X%*%beta[it,]
      cor <- as.vector(t(X)%*%res)
      if (max(abs(cor))>10^-8){
        S <- sign(cor)
        max.cor <- max(abs(cor))
        obj <- get.h(X,Y,beta[it,],as.vector(cor),active,it,XX,NZ)
        NZ <- obj$NZ
        h <- obj$h
        active <- active & (1:p)!=obj$slack
        act <- min((1:p)[active])
        direc.1 <- t(h)%*%(XX[act,]-XX)
        direc.2 <- t(h)%*%(XX[act,]+XX)
        dist.hitzero <- -beta[it,]/h
        index2 <- (1:p)[h!=0][dist.hitzero[h!=0]>10^-11]
        dist1 <- (cor[act]-cor)/direc.1
        dist2 <- (cor[act]+cor)/direc.2
        index <- c((1:p)[!active][dist1[!active]>10^-11],
        (1:p)[!active][dist2[!active]>10^-11])
        index3 <- c((1:p)[!active][dist1[!active]>10^-11], ((p+1):(2*p))[!active][dist2[! active]>10^-11])
        dist <- c(dist1,dist2)[index3]
        if (length(index)==0)
          gamma1 <- 1/0
        else
          gamma1 <- min(dist)
        if (length(index2)==0)
          gamma2 <- 1/0
        else
          gamma2 <- min(dist.hitzero[index2])
        gamma3 <- min(c(dist1[act],dist2[act]),na.rm=T)
        gamma <- min(gamma1,gamma2,gamma3)
        if (gamma1<gamma2)
          active[index[which.min(dist)]] <- T
        else{NNZ <- (1:p)[index2][which.min(dist.hitzero[index2])]
              NZ[c(NNZ,NNZ+p)] <- F}
        beta <- rbind(beta,beta[it,]+h*gamma)
        it <- it+1
      }
    }
    beta=t(t(beta)/normx)
    for (j in 1:p)
    beta[,j] <-beta[,j]/w[j]
    beta <-(1+lambda2)*(beta)
    obj <- list(beta=beta,mu=meany,meanx=meanx,normx=normx,type="LASSO")
    class(obj) <- "lars"
    obj
  }
