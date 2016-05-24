goric_penalty <-
function(object, iter=100000, type="GORIC", mc.cores=1){
  if (!(inherits(object, "orlm") | inherits(object, "orgls"))) stop("object needs to be of class orlm or orgls")
  if (all(object$constr == 0) & object$nec == 0 & type == "GORIC"){
    penalty <- length(object$coefficients) + 1
  } else {
    if (iter < 1) stop("No of iterations < 1")
    Sigma <- object$sigma
    x <- object$X
    invW <- kronecker(solve(Sigma), t(x) %*% x)
    W <- solve(invW)
    Z <- rmvnorm(n=iter, mean=rep(0, ncol(W)), sigma=W)
    Dmat <- 2*invW
    Amat <- t(object$constr)
    bvec <- object$rhs
    nec <- object$nec
    if (mc.cores == 1){
      nact <- apply(Z, 1, function(z){
        dvec <- 2*(z %*% invW)
        QP <- solve.QP(Dmat,dvec,Amat,bvec=bvec, meq=nec)
        if (QP$iact[1] == 0) return(0) else return(length(QP$iact))
      })
    } else {
      cl <- makeCluster(mc.cores)
      clusterExport(cl, c("solve.QP"))
      nact <- parRapply(cl, Z, function(z, invW,Dmat,Amat,bvec,nec){
        dvec <- 2*(z %*% invW)
        QP <- solve.QP(Dmat,dvec,Amat,bvec=bvec, meq=nec)
        if (QP$iact[1] == 0) return(0) else return(length(QP$iact))
      }, invW=invW, Dmat=Dmat, Amat=Amat, bvec=bvec, nec=nec)
      stopCluster(cl)
    }
    dimsol <- ncol(W) - nact
    LP <- sapply(1:(ncol(W)+1), function(x) sum(x == (dimsol+1)))/iter
    if (type == "GORIC"){
      penalty <- 1 + sum((1:ncol(W))*LP[-1])
    } else {
      N <- nrow(object$X)
      t <- ncol(object$y)
      k <- ncol(object$X)
      plp <- sum((1:ncol(W))*LP[-1])
      qlp <- -2*plp + sum(((1:ncol(W))^2)*LP[-1]) - plp^2
    }
    if (type == "GORICCa"){
      penalty <- -0.5*t*N + (0.5*t*N) * (t*N * ((t*N-plp)^2 + 2*t*N + qlp)) / ((t*N - plp)^3) + 0.5*sum((1:ncol(W))*LP[-1]*((t*N) / (t*N-(1:ncol(W)) - 2)))
    }
    if (type == "GORICCb"){
      penalty <- sum(((t*N * (0:ncol(W) + 1)) / (t*N - 0:ncol(W) - 2)) * LP)
    }
  }
  return(penalty)
}


orglm_penalty <- function(object, iter=100000, type="GORIC", mc.cores=1){
  if (!(inherits(object, "orglm"))) stop("object needs to be of class orglm")
  X <- object$X
  wt <- sqrt(object$weights)
  wX <- X * wt
  invW <- t(wX) %*% wX
  vc <- qr.solve(invW)
  Z <- rmvnorm(n=iter, mean=rep(0, ncol(vc)), sigma=vc)
  Dmat <- 2*invW
  Amat <- t(object$constr)
  bvec <- object$rhs
  nec <- object$nec
  if (mc.cores == 1){
    nact <- apply(Z, 1, function(z){
      dvec <- 2*(z %*% invW)
      QP <- solve.QP(Dmat,dvec,Amat,bvec=bvec, meq=nec)
      if (QP$iact[1] == 0) return(0) else return(length(QP$iact))
    })
  } else {
    cl <- makeCluster(mc.cores)
    clusterExport(cl, c("solve.QP"))
    nact <- parRapply(cl, Z, function(z, invW,Dmat,Amat,bvec,nec){
      dvec <- 2*(z %*% invW)
      QP <- solve.QP(Dmat,dvec,Amat,bvec=bvec, meq=nec)
      if (QP$iact[1] == 0) return(0) else return(length(QP$iact))
    }, invW=invW, Dmat=Dmat, Amat=Amat, bvec=bvec, nec=nec)
    stopCluster(cl)
  }
  dimsol <- ncol(vc) - nact
  LP <- sapply(1:(ncol(vc)+1), function(x) sum(x == (dimsol+1)))/iter
  if (type == "GORIC"){
    penalty <- sum((1:ncol(vc))*LP[-1])
  } else {
    N <- nrow(object$X)    
    k <- ncol(object$X)
    plp <- sum((1:ncol(vc))*LP[-1])
    qlp <- -2*plp + sum(((1:ncol(vc))^2)*LP[-1]) - plp^2
  }
  if (type == "GORICCa"){
    penalty <- -0.5*N + (0.5*N) * (N * ((N-plp)^2 + 2*N + qlp)) / ((N - plp)^3) + 0.5*sum((1:ncol(vc))*LP[-1]*((N) / (N-(1:ncol(vc)) - 2)))
  }
  if (type == "GORICCb"){
    penalty <- sum(((N * (0:ncol(vc) + 1)) / (N - 0:ncol(vc) - 2)) * LP)
  }
  return(penalty)
}


