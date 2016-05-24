dNNiW <- function(xi, psi, Sigma, B, hyperprior, log=FALSE){
  
  res=NA
  
  xi0 <- hyperprior[["b_xi"]]
  psi0 <- hyperprior[["b_psi"]]
  nu0 <- hyperprior[["nu"]]
  B0 <- hyperprior[["B"]]
  if(is.null(B0)){
    D0_xi <- hyperprior[["D_xi"]]
    D0_psi <- hyperprior[["D_psi"]]
    B0 <- diag(c(D0_xi, D0_psi))
  }else{
    D0_xi <- B0[1,1]
    D0_psi <- B0[2,2]
  }
  lambda0 <- hyperprior[["lambda"]]
  
  p <- ncol(lambda0)
  
  if(!is.list(Sigma)){
    if(length(dim(Sigma))==3){
      Sigma <- apply(X=Sigma, MARGIN=3, FUN=list)
      Sigma <- lapply(Sigma, FUN='[[', 1)
      B <- apply(X=B, MARGIN=3, FUN=list)
      B <- lapply(B, FUN='[[', 1)
      xi <- apply(X=xi, MARGIN=2, FUN=list)
      xi <- lapply(xi, FUN='[[', 1)
      psi <- apply(X=psi, MARGIN=2, FUN=list)
      psi <- lapply(psi, FUN='[[', 1)
    }
    else{
      B <- list(Sigma)
      Sigma <- list(Sigma)
      xi <- list(xi)
      psi <- list(psi)            
    }
  }
  
  mu0 <- c(xi0, psi0)
  D0 <- c(rep(D0_xi, length(xi0)),rep(D0_psi, length(psi0)))
  mu <- mapply(function(a, b){c(a,b)}, a=xi, b=psi, SIMPLIFY=F)
  S0 <- B0%x%lambda0
  S <- mapply('%x%',X=B, Y=Sigma, SIMPLIFY=F)
  
  try_inv_<-mapply(function(S) {try(solve(S),silent=TRUE)},S=S,SIMPLIFY=F)
  try_inv<-which(mapply(function(q){inherits(q,"try-error")},q=try_inv_))
  if (length(try_inv)){
    for (elem in try_inv){
      S[[elem]]<-  S[[elem]]+diag(ncol(S[[1]])) 
    }
  }
  if(!log){
    res <- mapply(function(mu, S){
      det(S)^(-(nu0+p+2)/2)*
        exp(-sum(diag(S0%*%solve(S)))/2
            -((mu-mu0)/D0)%*%crossprod(solve(S), (mu-mu0))/2)
    }, mu=mu, S=S)
  }else{
    res <- mapply(function(mu, S){
      -(nu0+p+2)/2*log(det(S))
      -sum(diag(S0%*%solve(S)))/2
      -((mu-mu0)/D0)%*%crossprod(solve(S), (mu-mu0))/2
    }, mu=mu, S=S)
  }
  
  return(res)
}