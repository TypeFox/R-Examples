
uniprobit <- function(mu,dmu,S,dS,y,w=NULL,indiv=FALSE,...) {
  sigma <- S^0.5
  alpha <- alpha0 <- pnorm(mu,sd=sigma)
  alpha[y==0] <- 1-alpha[y==0]
  M <- -sigma*dnorm(mu/sigma)
  V <- S*alpha0+mu*M  
  U1 <- 0.5*t(dS)%*%(-alpha0+V/S)/S
  U <- matrix(0,ncol(dmu)+nrow(U1),ncol(U1))
  if (is.null(w)) w <- rep(1,ncol(U))
  for (i in seq(ncol(U))) {
    U[,i] <- c(- dmu[i,,drop=FALSE]*(1/S*M[i]),U1[,i])/alpha[i]*
       ifelse(y[i]==0,-1,1)*w[i]
  }  
  if (indiv) return(structure(t(U),logLik=log(alpha)))
  return(structure(rowSums(U),logLik=sum(log(alpha))))
}
