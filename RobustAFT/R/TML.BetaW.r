TML.BetaW <- function(X,y,delta,Beta,sigma,Beta.t,sigma.t,cl,cu,maxit,tol,nitmon)
{
# iteratively reweighting algorithm for Beta in TML log Weibull
p <- length(Beta); n <- length(y); nu <- sum(delta); nc <- n-nu; zero <- 1e-6
nit <- 1; Beta1 <- rep(100,p)
indu <- (1:n)[delta==1]; indc <- (1:n)[delta==0]
mui.t <- X %*% as.matrix(Beta.t)
rs.t  <- (y-mui.t)/sigma.t
wgt   <- ww(rs.t, cl, cu)
#
while ( max(abs(Beta1-Beta)) > tol & (nit < maxit) ) {
D1  <- D2 <- ym <- rep(0,n); vi <- ti <- rep(0,nc)
nit <- nit+1; Beta1 <- Beta
mui <- X %*% as.matrix(Beta1)
rs  <- (y-mui)/sigma
gi  <- rep(1,nu)
if (nu > 0) {
  ru       <- rs[indu]
  cnd      <- ru != 0
  gi[cnd]  <- ps0(ru[cnd])/ru[cnd]
  D1[indu] <- wgt[indu]*gi
  D2[indu] <- wgt[indu]*gi
  ym[indu] <- y[indu]}
if (nu < n) {
  I1 <- I0 <- rep(0,nc)
  rc <-  rs[delta==0]
  muic <- mui[indc]
  muit <- mui.t[indc]
  ai <- pmax( rc, (sigma.t*cl - muic + muit )/sigma   )
  bi <-           (sigma.t*cu - muic + muit )/sigma 
  den <- 1-plweibul(rc)
  ok  <- den > zero    
  for (i in 1:nc) { 
    if (bi[i] > ai[i]) I0[i] <- integrate(intg0.TMLW,ai[i],bi[i])$val    
    if (bi[i] > ai[i]) I1[i] <- dlweibul(ai[i])-dlweibul(bi[i])  }
  I1 <- sigma*I1+muic*I0    
  vi[ok] <- I0[ok]/den[ok]
  ti[ok] <- I1[ok]/den[ok] 
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