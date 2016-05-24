TML.BetaG <- function(X,y,delta,Beta,sigma,Beta.t,sigma.t,cu,maxit,tol,nitmon)
{
  cl <- -cu
# iteratively reweighting algorithm for Beta in TML Gauss 
p <- length(Beta); n <- length(y); nu <- sum(delta); nc <- n-nu
nit <- 1; Beta1 <- rep(100,p); zero <- 1e-6
indu <- (1:n)[delta==1]; indc <- (1:n)[delta==0]
mui.t <- X %*% as.matrix(Beta.t)
rs.t  <- (y-mui.t)/sigma.t
wgt   <- ww(rs.t, cl, cu)
while ( max(abs(Beta1-Beta)) > tol & (nit < maxit) ) {
D1  <- D2 <- ym <- rep(0,n); vi <- ti <- rep(0,nc)
nit <- nit+1; Beta1 <- Beta
mui <- X %*% as.matrix(Beta1)
rs  <- (y-mui)/sigma
if (nu > 0) {
  D1[indu] <- wgt[indu]
  D2[indu] <- wgt[indu]
  ym[indu] <- y[indu]}
if (nu < n) {
  I1 <- I0 <- rep(0,nc)
  rc <-  rs[delta==0]
  muic <- mui[indc]
  muit <- mui.t[indc]
  ai <- pmax( rc, (sigma.t*cl - muic + muit )/sigma   )
  bi <-             (sigma.t*cu - muic + muit )/sigma
  den <- 1-pnorm(rc)
  ok  <- den > zero
   for (i in 1:nc) {     
    if (bi[i] > ai[i]) I0[i] <- pnorm(bi[i])-pnorm(ai[i])
    if (bi[i] > ai[i]) I1[i] <- dnorm(ai[i])-dnorm(bi[i]) }
  I1 <- sigma*I1+muic*I0
  vi[ok] <- I0[ok]/den[ok]
  ti[ok] <- I1[ok]/den[ok] 
  D1[indc] <- vi
  D2[indc] <- ti
  ym[indc] <- rep(1,nc)}
A <- t(X)%*%(D2*ym)
B <- t(X)%*%(D1*X)
Beta <- solve(B)%*%A
if(nitmon) cat(nit, Beta, Beta1, "\n")}
list(Beta=Beta,nit=nit)}

