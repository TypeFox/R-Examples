dNiW <- function(mu, Sigma, hyperprior, log=FALSE){
    
    res=NA
    
    mu0 <- hyperprior[["mu"]]
    nu0 <- hyperprior[["nu"]]
    kappa0 <- hyperprior[["kappa"]]
    lambda0 <- hyperprior[["lambda"]]
    
    p <- ncol(lambda0)
    
    if(!is.list(Sigma)){
        Sigma <- apply(X=Sigma, MARGIN=3, FUN=list)
        Sigma <- lapply(Sigma, FUN='[[', 1)
        mu <- apply(X=mu, MARGIN=2, FUN=list)
        mu <- lapply(mu, FUN='[[', 1)
    }
    
    if(!log){
    res <- mapply(function(mu, Sigma){
        det(Sigma)^(-(nu0+p+2)/2)*
            exp(-sum(diag(lambda0%*%solve(Sigma)))/2
                -kappa0/2*(mu-mu0)%*%crossprod(solve(Sigma), (mu-mu0)))
        }, mu=mu, Sigma=Sigma)
    }else{
      res <- mapply(function(mu, Sigma){
        -(nu0+p+2)/2*log(det(Sigma))
        -sum(diag(lambda0%*%solve(Sigma)))/2
        -kappa0/2*(mu-mu0)%*%crossprod(solve(Sigma), (mu-mu0))
      }, mu=mu, Sigma=Sigma)
    }
    
    return(res)
}