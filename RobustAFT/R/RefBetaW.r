RefBetaW <- function(X,y,delta,Beta,sigma,maxit,tol,nitmon)
{
# iteratively rewighting algorithm for Beta refinement, log-Weibull case 
p <- length(Beta); n <- length(y); nu <- sum(delta); nc <- n-nu; c0 <- -0.1351788
nit <- 1; Beta1 <- rep(100,p); zero <- 1e-6
indu <- (1:n)[delta==1]
indc <- (1:n)[delta==0]
#
while ( max(abs(Beta1-Beta)) > tol & (nit < maxit) ) {
I1 <- I0 <- vi <- ti <- rep(0,nc)
D1 <- D2 <- ym <- rep(0,n)
nit <- nit+1; Beta1 <- Beta
mui <- X %*% as.matrix(Beta1)
rs  <- (y-mui)/sigma
if (nu > 0) {
  ru      <- rs[delta==1]
  wi      <- PspSw(ru-c0)
  cnd     <- ru != 0
  wi[cnd] <- PsiSw(ru[cnd]-c0)/(ru[cnd]-c0)  }
if (nu < n) {
  rc <-  rs[delta==0]
  den <- 1-plweibul(rc)
  ok  <- den > zero
  for (i in 1:nc) {
     I0[i] <- integrate(intw0,lower=rc[i]-c0,upper=1.718)$val
     I1[i] <- integrate(intw1,lower=rc[i]-c0,upper=1.718)$val}
  I1 <- sigma*I1+mui[delta==0]*I0    
  vi[ok] <- I0[ok]/den[ok]
  ti[ok] <- I1[ok]/den[ok] }
if (nu > 0) {
  D1[indu] <- wi
  D2[indu] <- wi
  ym[indu] <- y[indu]-sigma*c0}
if (nu < n) {
  D1[indc] <- vi
  D2[indc] <- ti
  ym[indc] <- rep(1,nc)}
A <- t(X)%*%(D2*ym)
B <- t(X)%*%(D1*X)
Beta <- solve(B)%*%A
    if(nitmon) cat(nit, Beta, Beta1, "\n")
  }
  list(Beta=Beta,nit=nit)
}