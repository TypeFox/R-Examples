pcovr_est <-
  function(X,Y,r,a,cross=FALSE,fold="LeaveOneOut"){
    N <- nrow(X)
    if (is.vector(Y)){
      if (length(Y) != N){
        print('The number of observations is not identical for X and Y')
      }
    } else if (N!=nrow(Y)){
      print('The number of observations is not identical for X and Y')
    }
    
    #compute projector on X
    S <- t(X) %*% X
    if (det(S)<1e-12){ 
      l <- eigen(S)$values
      w <- which(l>1e-8*max(l))
      Sh <- eigen(S)$vectors[,w] %*% diag(l[w]^-1) %*% t(eigen(S)$vectors[,w])
      Hx <- X %*% Sh %*% t(X)
    } else {
      Hx <- X %*% solve(S) %*% t(X)
    }
    
    # Compute PCovR solution
    block1 <- sqrt(a/SUM(X)$ssq)*X
    block2 <- Hx %*% Y
    block2 <- sqrt((1-a)/SUM(Y)$ssq)*block2
    V <- cbind(block1,block2)
    P <- V %*% eigen(t(V) %*% V)$vectors %*% diag(eigen(t(V) %*% V)$values^(-1/2))
    Te <- P[,1:r]
    W <- ginv(X) %*% Te
    Px <- t(Te) %*% X
    Py <- t(Te) %*% Y
    Rx2 <- SUM(Te %*% Px)$ssq*SUM(X)$ssq^-1
    Ry2 <- SUM(Te %*% Py)$ssq*SUM(Y)$ssq^-1
    B <- W %*% Py
    
    if (cross==TRUE){
      if (fold == "LeaveOneOut"){
        fold <- N
      }
      Yhatcv <- Y
      count <- 1
      LeaveOut <- rep(floor(N/fold), times = fold)
      LeaveOut[1] <- LeaveOut[1] + N%%fold
      CumulSum <- cumsum(LeaveOut)
      for (i in 1:fold){
        Xi <- X
        Xi <- Xi[-(count:CumulSum[i]),]
        Yi <- Y
        Yi <- array(Yi[-(count:CumulSum[i]),],c(nrow(X)-LeaveOut[i],ncol(Y)))
        
        S <- t(Xi) %*% Xi
        if (det(S)<1e-12){ 
          l <- eigen(S)$values
          w <- which(l>1e-8*max(l))
          Sh <- eigen(S)$vectors[,w] %*% diag(l[w]^-1) %*% t(eigen(S)$vectors[,w])
          Hx <- Xi %*% Sh %*% t(Xi)
        } else {
          Hx <- Xi %*% solve(S) %*% t(Xi)
        }
        
        block1 <- sqrt(a/SUM(Xi)$ssq)*Xi
        block2 <- Hx %*% Yi
        block2 <- sqrt((1-a)/SUM(Yi)$ssq)*block2
        V <- cbind(block1,block2)
        P <- V %*% eigen(t(V) %*% V)$vectors %*% diag(eigen(t(V) %*% V)$values^(-1/2))
        Ti <- P[,1:r]
        Wi <- ginv(Xi) %*% Ti
        Bi <- Wi %*% t(Ti) %*% Yi
        Yhatcv[(count:CumulSum[i]),] <- X[(count:CumulSum[i]),] %*% Bi
        count <- count + LeaveOut[i]
      }
      Qy2 <- 1-SUM(Y-Yhatcv)$ssq*SUM(Y)$ssq^-1
      results <- list(W=W,B=B,Rx2=Rx2,Ry2=Ry2,Te=Te,Px=Px,Py=Py,Qy2=Qy2)
    } else {
      results <- list(W=W,B=B,Rx2=Rx2,Ry2=Ry2,Te=Te,Px=Px,Py=Py)
    }
    return(results)
  }