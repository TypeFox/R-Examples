predict.pfda <-
function(object,X,...){
  # Initialization
  Y = as.matrix(X)
  n = nrow(Y)
  p = ncol(Y)
  K = object$prms$K
  U = object$V
  mu = object$prms$mean
  prop = object$prms$prop
  D = object$prms$D
  d = ncol(U)
  # browser()
  QQ = matrix(NA,n,K)
  QQ2 = matrix(NA,n,K)
  T = matrix(NA,n,K)
  
  # Compute posterior probabilities
  for (k in 1:K){
    
    bk = D[k,p,p]
    mY = object$prms$my[k,]
    YY = Y-t(matrix(rep(mY,n),ncol(Y),nrow(Y)))
    projYY = YY %*% U %*% t(U) # proj dans l'espace latent
    
    if (d==1){
      for (i in 1:n){QQ[i,k] =  1/D[k,1,1] * sum(projYY[i,]^2) + 1/D[k,p,p]*sum((YY[i,] - projYY[i,])^2) + (p-d)*log(bk) + log(D[k,1,1]) - 2*log(prop[k])}
    }
    else{
      sY = U %*% ginv(D[k,(1:d),(1:d)]) %*% t(U)
      for (i in 1:n){QQ[i,k] =  projYY[i,] %*% sY %*% as.matrix(projYY[i, ],p,1) + 1/bk*sum((YY[i,] - projYY[i,])^2) + (p-d)*log(bk) + log(det(D[k,(1:d),(1:d)])) - 2*log(prop[k])}
    }
  }
  
  for (k in 1:K) {T[,k] = 1 / rowSums(exp(0.5*(t(.pfda.repmat(QQ[,k],K,1))-QQ)))}
  list(P=T,cls=max.col(T))
}
