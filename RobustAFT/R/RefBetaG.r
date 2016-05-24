RefBetaG <- function(X,y,delta,Beta,sigma,maxit,tol,nitmon) {
zero <- 1e-6
# iteratively rewighting algorithm for Beta refinement, Gaussian case 
p <- length(Beta); n <- length(y); nu <- sum(delta); nc <- n-nu
nit <- 1; Beta1 <- rep(100,p)
indu <- (1:n)[delta==1]
indc <- (1:n)[delta==0]
while ( max(abs(Beta1-Beta)) > tol & (nit < maxit) ) {
I1 <- I0 <- I2 <- vi <- ui <- wi <- rep(0,nc)
D1 <- D2 <- ym <- rep(0,n)
nit <- nit+1; Beta1 <- Beta
mui <- X %*% as.matrix(Beta1)
rs  <- (y-mui)/sigma
if (nu > 0) {
  ru      <- rs[delta==1]
  wi      <- PspSG(ru)
  cnd     <- ru != 0
  wi[cnd] <- (PsiSG(ru[cnd])/ru[cnd])  }
if (nu < n) {
  rc  <-  rs[delta==0]
  den <- 1-pnorm(rc)
  ok  <- den > zero
  for (i in 1:nc) {
     rci <- rc[i]
     if (rci < -6) rci <- -6
     I0[i] <- integrate(intg0,lower=rci,upper=6)$val
#    cat("i, rc_i, I0_i, I2_i:",i,rc[i],I0[i],I2[i],"\n")
     I2[i] <- integrate(intg2,lower=rci,upper=6)$val}
  I1 <- sigma*I2+mui[delta==0]*I0
  vi[ok] <- I1[ok]/den[ok]
  ui[ok] <- I0[ok]/den[ok] }
D1[indu] <- wi
D1[indc] <- vi
D2[indu] <- wi
D2[indc] <- ui
ym[indu] <- y[indu]
ym[indc] <- rep(1,nc)
A <- t(X)%*%(D1*ym)
B <- t(X)%*%(D2*X)
Beta <- solve(B)%*%A
if(nitmon) cat(nit, Beta, Beta1, "\n")
  }
  list(Beta=Beta,nit=nit)
}
