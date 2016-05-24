bwchoice <- function(X,objectif,kernelx="g",itermax=1000) {
  p <- ncol(X)
  res <- rep(0,p)
  if (length(objectif)==1) objectif <- rep(objectif,p)
  poids1d <- function(kernelx="g",X,bx,valx){
    vecteur <- (X-valx)/bx
    if(kernelx=="e"){noyau <- epane}
    if(kernelx=="g"){noyau <- gaussien}
    if(kernelx=="q"){noyau <- quartic}
    if(kernelx=="u"){noyau <- uniform}
    vv <- noyau(vecteur)
    return(W=vv)
  }
  H1d0 <- function(kernelx="g",X,bx){
    n <- length(X)
    H <- matrix(0,ncol=n,nrow=n)
    for(i in 1:n){
      w <- poids1d(kernelx,X,bx,X[i])
      H[i,] <- w/sum(w)
    }
    trace <- sum(diag(H))
    return(trace)}

  choixddlparvar <- function(fenetre,X,objectif,kernelx) {
    prov <- H1d0(X=X,bx=fenetre,kernelx =kernelx)
    res <- prov-objectif
    return(res)
  }
  if (any(objectif<=1)) stop("degree of freedom should be greater than 1\n")
  for (j in 1:ncol(X)) {
    depart <- 3*abs(diff(range(X[,j])))
    if (choixddlparvar(depart,X[,j],objectif[j],kernelx)>0) {
      while (choixddlparvar(depart,X[,j],objectif[j],kernelx)>0) {
        depart <- depart*2
      }
    }
    res[j] <- uniroot(choixddlparvar,interval=c(depart,1e-10),X=X[,j],objectif=objectif[j],maxiter=itermax,kernelx=kernelx)$root
  }
  return(res)
}
