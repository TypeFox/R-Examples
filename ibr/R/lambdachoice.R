lambdachoice <- function(X,ddlobjectif,m=2,s=0,itermax,smoother="tps") {
  n <- nrow(X)
  d <- ncol(X)
  p <- 2*m-d
  ddlmin <- choose(m+d-1,m-1)
  if (ddlobjectif<=ddlmin) stop(paste("the objective df is too small, choose it greater than",ddlmin))
  Sgu <- DuchonS(X,m)
  qrSgu <- qr(Sgu)
  F2 <- qr.Q(qrSgu,complete=TRUE)[,-(1:ncol(Sgu))]
  Kgu <- DuchonQ(X,0,m,s,symmetric=TRUE)
  ainv <- t(F2)%*%Kgu%*%F2
  vp <- eigen(ainv,only.values=TRUE,symmetric=TRUE)$values
  trace <- function(loglambda,vp1=vp) {
    n-sum(1/(vp1/exp(loglambda)+1)) - ddlobjectif
  }
  l1 <- 1
 for (k in 1:25) {
    tr <- n-sum(1/(1 + vp/l1))
       if (tr <= ddlobjectif)  break
        l1 <- l1 * 4
    }
     l2 <- 1
    for (k in 1:25) {
        tr <- n-sum(1/(1 + vp/l2))
         if (tr > ddlobjectif) break
        l2 <- l2/4
    }
    resultat <- uniroot(trace,c(log(l2),log(l1)),vp1=vp,maxiter =itermax)
    return(exp(resultat$root))
}
