HWChisqStats <- function(X,x.linked=FALSE,pvalues=FALSE) {
  if(!x.linked) {
    n <- rowSums(X)
    nA <- 2*X[,1] + X[,2]
    nB <- 2*X[,3] + X[,2]
    chi <- (4*X[,1]*X[,3] - X[,2]*X[,2])/(nA*nB)
    if(pvalues) stat <- pchisq(n*chi*chi,1,lower.tail=FALSE) else stat <- n*chi*chi
  } else {
    nm <- X[,1] + X[,2]
    nA <- X[,1] + 2*X[,3] + X[,4]
    nA2 <- nA*nA
    nf <- X[,3] + X[,4] + X[,5]
    nt <- nm + 2*nf
    nB <- nt - nA
    nB2 <- nB*nB
    m.num <- nA*X[,2]-nB*X[,1]
    X.male <- (m.num*m.num)/(nm*nA*nB)
    f.num <- nt*nt*(X[,3]*X[,3]*nB2 + 0.5*X[,4]*X[,4]*nA*nB + X[,5]*X[,5]*nA2) - nf*nf*nA2*nB2
    X.female <- f.num/(nf*nA2*nB2)
    if(pvalues) stat <- pchisq(X.male+X.female,2,lower.tail=FALSE) else stat <- X.male+X.female
  }
}
