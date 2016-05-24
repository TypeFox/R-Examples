##
## alrtest
##
alrtest <- function(z, A, r){
  if(!(class(z)=="ca.jo")){
    stop("\nPlease, provide object of class 'ca.jo' as 'z'.\n")
  }
  r <- as.integer(r)
  A <- as.matrix(A)
  if(!(nrow(A)==z@P)){
    stop("\nRow number of 'A' is unequal to VAR order.\n")
  }
  if(r >= z@P || r<1){
    stop("\nCount of cointegrating relationships is out of allowable range.\n")
  }
  type <- "Estimation and testing under linear restrictions on beta"
  B <- qr.Q(qr(A), complete=TRUE)[,-c(1:ncol(A))]
  N <- nrow(z@Z0)
  M00 <- crossprod(z@Z0)/N
  M11 <- crossprod(z@Z1)/N
  MKK <- crossprod(z@ZK)/N
  M01 <- crossprod(z@Z0, z@Z1)/N
  M0K <- crossprod(z@Z0, z@ZK)/N
  MK0 <- crossprod(z@ZK, z@Z0)/N
  M10 <- crossprod(z@Z1, z@Z0)/N
  M1K <- crossprod(z@Z1, z@ZK)/N
  MK1 <- crossprod(z@ZK, z@Z1)/N
  M11inv <- solve(M11)
  S00 <- M00 - M01%*%M11inv%*%M10
  S0K <- M0K - M01%*%M11inv%*%M1K
  SK0 <- MK0 - MK1%*%M11inv%*%M10
  SKK <- MKK - MK1%*%M11inv%*%M1K
  Sab <- t(A)%*%S00%*%B
  Skb <- t(S0K)%*%B
  Sbb <- t(B)%*%S00%*%B
  Sbbinv <- solve(Sbb)
  RA <- z@R0%*%A - z@R0%*%B%*%Sbbinv%*%t(Sab)
  RK <- z@RK - z@R0%*%B%*%Sbbinv%*%t(Skb)
  Saa.b <- crossprod(RA, RA)/N
  Sak.b <- crossprod(RA, RK)/N
  Ska.b <- crossprod(RK, RA)/N
  Skk.b <- crossprod(RK, RK)/N
  Ctemp <- chol(Skk.b, pivot=TRUE)
  pivot <- attr(Ctemp, "pivot")
  oo <- order(pivot)
  C <- t(Ctemp[,oo])
  Cinv <- solve(C)
  Saa.binv <- solve(Saa.b)
  valeigen <- eigen(Cinv%*%Ska.b%*%Saa.binv%*%Sak.b%*%t(Cinv), symmetric = TRUE)
  lambda.res <- valeigen$values
  e <- valeigen$vector
  V <- t(Cinv)%*%e
  V <- as.matrix(V[,1:r])
  Vorg <- V
  idx <- 1:r
  V <- sapply(idx, function(x) V[ , x] / V[1,x])
  PHI <- solve(t(A)%*%A)%*%Sak.b%*%Vorg
  ALPHA <- as.matrix(A%*%PHI)
  ALPHAorg <- ALPHA
  ALPHA <- sapply(idx, function(x) ALPHA[ , x] * Vorg[1,x])
  PI <- ALPHA %*% t(V)
  GAMMA <- M01%*%M11inv - PI%*%MK1%*%M11inv
  DELTA.bb <- Sbb
  DELTA.ab <- Sab - t(A)%*%ALPHA%*%t(V)%*%Skb
  DELTA.aa.b <- Saa.b - t(A)%*%ALPHA%*%t(ALPHA)%*%A
  lambda <- z@lambda
  teststat <- N*sum(log((1-lambda.res[1:r])/(1-lambda[1:r])))
  df <- r*(z@P - ncol(A))
  pval <- c(1-pchisq(q = teststat, df = df), df)
  new("cajo.test", Z0=z@Z0, Z1=z@Z1, ZK=z@ZK, ecdet=z@ecdet, H=NULL, A=A, B=B, type=type, teststat=teststat, pval=pval, lambda=lambda.res, Vorg=Vorg, V=V, W=ALPHA, PI=PI, DELTA=NULL, DELTA.bb=DELTA.bb, DELTA.ab=DELTA.ab, DELTA.aa.b=DELTA.aa.b, GAMMA=GAMMA, test.name="Johansen-Procedure")
}
