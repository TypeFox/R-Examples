# ML_functions.s

# Auxiliary functions for (weighted) maximum likelihood
# ---------------------------------------------------

sigma.mleqn.loggamma <- function(sigma,mu,lam,d,yi,wi){
# LHS of ML equation for sigma
n  <- length(yi); swi <- sum(wi)
rs <- (yi-mu)/sigma
sum(wi*csiLG(rs,lam)*rs)/swi - (d-1)}

lam.mleqn.loggamma <- function(lambda,mu,sigma,aa,yi,wi){
# LHS of ML equation for lambda
n  <- length(yi); swi <- sum(wi)
rs <- (yi-mu)/sigma
sum(wi*psiLG(rs,lambda) )/swi - aa}

Fxdsig.loggamma <- function(sig0,mu0,lam0,yi,wi,d,tol=0.001,maxit=100){
# Fixed point step for sigma
sig <- sig0; nit <- 1
repeat{
  sig0 <- sig
  rs   <- (yi-mu0)/sig
  sig  <- sum(wi*csiLG(rs,lam0)*rs)/sum(wi)*sig/(d-1)
  delta <- sig-sig0
  if (nit==maxit | abs(delta) < tol) break
  nit <- nit+1}
list(sig=sig,nit=nit,delta=delta)}

Nwtlam.loggamma <- function(sig0,mu0,lam0,yi,wi,a,lam.low,lam.sup,tol=0.001,maxit=50){
# Newton step for lambda
lam   <- lam0; nit <- 1; swi <- sum(wi)
rs    <- (yi-mu0)/sig0
repeat {
 f    <- lam.mleqn.loggamma(lam,mu0,sig0,a,yi,wi)
 f1   <- sum(wi*psiLG.dot(rs,lam) )/swi
 if (!is.finite(f) | !is.finite(f1))    {delta <- NA; lam <- NA; break}
 delta  <- -f/f1
 if (is.nan(delta) | !is.finite(delta)) {delta <- NA; lam <- NA; break}
 lam <- lam+delta
 if (lam < lam.low | lam > lam.sup )    {delta <- NA; lam <- NA; break}
 if (nit==maxit | abs(delta) < tol) break
 nit <- nit+1}
list(lam=lam,nit=nit,delta=delta)}

Nwtsig.loggamma  <- function(sig0,mu0,lam0,yi,wi,d,tol=0.001,maxit=50){
# Fixed point step for lambda
sig <- sig0; nit <- 1; swi <- sum(wi)
repeat{ 
  sig0 <- sig
  rs <- (yi-mu0)/sig
  f  <- sigma.mleqn.loggamma(sig,mu0,lam0,d,yi,wi)
  f1 <- sum(wi*rs*(dersig.csiLG(rs,lam0,sig)-csiLG(rs,lam0)/sig))/swi
  if (!is.finite(f) | !is.finite(f1))    {delta <- NA; sig <- NA; break}
  delta  <- -f/f1
  if (is.nan(delta) | !is.finite(delta)) {delta <- NA; sig <- NA; break}
  sig <- sig+delta
  if (sig<=0) {delta <- NA; sig <- NA; break}
  if (nit==maxit | abs(delta) < tol) break
  nit <- nit+1}
list(sig=sig,nit=nit,delta=delta)}

Rgfsig.loggamma <- function(sig0,mu0,lam0,yi,wi,dd,tol=0.001,maxit=50,sig.low,sig.sup){
# Regula falsi for ML equation for sigma
sigl <- sig.sup
fun  <- sigma.mleqn.loggamma(sigl,mu0,lam0,dd,yi,wi)
if (is.na(fun)) fun <- -1000
if (fun < 0) {cat("Rgfsig: sigma.mleqn.loggamma(sig.sup) < 0","\n"); return(list(sig=-sig0,nit=0))}
while(fun > 0) {
  sigl     <- sigl/2
  if (sigl < sig.low) {cat("Rgfsig: sigl too small", "\n"); return(list(sig=-sig0,nit=0))}
  fun <- sigma.mleqn.loggamma(sigl,mu0,lam0,dd,yi,wi) }
sigu <- sigl
fun  <- sigma.mleqn.loggamma(sigu,mu0,lam0,dd,yi,wi) 
while(fun < 0) {
  sigu <- sigu+0.5
  if (sigu > sig.sup) {cat("Rgfsig: sigu too large", "\n"); return(list(sig=-sig0,nit=0))}
  fun <- sigma.mleqn.loggamma(sigu,mu0,lam0,dd,yi,wi) }
z    <- regfal(sigma.mleqn.loggamma,0,lower=sigl,upper=sigu,nint=10,tol=0.001,maxit=50,
               mu=mu0,lam=lam0,d=dd,yi=yi,wi=wi)
sig  <- z$solution; nit <- z$nit
list(sig=sig,nit=nit)}

Rgflam.loggamma <- function(sig0,mu0,lam0,yi,wi,aa,tol=0.001,maxit=50,lam.low,lam.sup){
# Regula falsi for ML equation for lambda
fun0 <- lam.mleqn.loggamma(lam0,   mu0,sig0,aa,yi,wi)
if (abs(fun0) < tol) return(list(lam=lam0,nit=0))
funl <- lam.mleqn.loggamma(lam.low,mu0,sig0,aa,yi,wi)
funu <- lam.mleqn.loggamma(lam.sup,mu0,sig0,aa,yi,wi)
laml <- lam.low; lamu <- lam.sup
if (fun0 > 0) {lamu <- lam0; funu <- fun0} else {laml <- lam0; funl <- fun0}
if (funl * funu > 0) {cat("Rgflam: no change of sign","\n"); return(list(lam=-lam0,nit=NA))}
# cat(laml,funl,lamu,funu,"\n")
z   <- regfal(lam.mleqn.loggamma,0,lower=laml,upper=lamu,nint=10,tol=tol,maxit=maxit,
              mu=mu0,sigma=sig0,a=aa,yi=yi,wi=wi)
lam <- z$solution; nit <- z$nit
list(lam=lam,nit=nit)}

# Weighted Maximum Likelihood (with fixed weights)

WML.loggamma <- function(y,wi,mu0,sig0,lam0,aa=0,dd=0,lam.low,lam.sup,tol=0.001,maxit=100,lstep=TRUE,sigstep=TRUE){
# solves ML equations for lambda in (lam.low,lam.sup); mu0, sig0,lam0 are the initial values; the weights wi are fixed
sig.low <- 1e-3*sig0; sig.sup <- 10000; nit <- 1
dsig <- dlam <- dmu <- 0; p <- 0; n <- length(y)
sig <- sig0; lam <- lam0; mu <- mu0
repeat {
# 
# scale step
  sig1 <- sig
  if (sigstep) {
    z <- Nwtsig.loggamma(sig1,mu,lam,y,wi,dd,tol,maxit=maxit)$sig
    if (!is.na(z)) sig <- z else {
      z   <- Rgfsig.loggamma(sig1,mu,lam,y,wi,dd,tol=0.001,maxit=50,sig.low=sig.low,sig.sup=sig.sup)
      sig <- z$sig
      if (sig == -sig1) return(list(mu=mu0,sig=sig0,lam=lam0,nit=NA)) }}
  dsig <- sig-sig1
# location step
  mu1     <- mu
  rs      <- (y-mu1)/sig
  di      <- rep(1,n)
  cnd     <- rs!=0
  di[cnd] <- -csiLG(rs[cnd],lam)/rs[cnd]
  mu      <- mean( wi*di*y )/mean( wi*di )
  dmu     <- mu - mu1
# shape step
  lam1 <- lam
  if (lstep) {
    z <- Nwtlam.loggamma(sig,mu,lam1,y,wi,aa,lam.low,lam.sup,tol,maxit=maxit)$lam
    if (!is.na(z))  lam <- z else {
      z <- Rgflam.loggamma(sig,mu,lam1,y,wi,aa,tol=0.001,maxit=50,lam.low,lam.sup)
      lam <- z$lam      
      if (lam == -lam1) return(list(mu=mu0,sig=sig0,lam=lam0,nit=NA))}}
  dlam <- lam-lam1
# cat(nit,lam,mu,sig,"\n")
  if (nit==maxit)  cat("WML: nit=maxit","\n")                                                
  if (nit==maxit | abs(dmu) < tol  &  abs(dsig) < tol & abs(dlam) < tol ) break
  nit <- nit+1}
list(mu=mu,sig=sig,lam=lam,nit=nit)}

# General regula falsi for solving f(x)=cc between lower and upper
# ----------------------------------------------------------------
regfal <- function(f,cc,lower,upper,nint=50,tol=0.001,maxit=20,...) {
  as <- seq(from=lower,to=upper,length=nint)
  ck <- NULL; nit <- 0; i <- 0
  while(i < (nint-1)) {
    i <- i+1
    fa <- f(as[i]  ,...)-cc
    fb <- f(as[i+1],...)-cc
# cat(i,as[i],as[i+1],fa,fb,"\n")
    if (fa*fb <= 0)
      break
  }
  if (i==(nint-1) & fa*fb>0) {
    cat("no solution ", "i=", i, "\n"); return(list(solution=ck,nit=-1))
  }
  ak <- as[i]
  bk <- as[i+1]
  fa <- f(ak,...)-cc
  fb <- f(bk,...)-cc
  nit <- 0; conv <- F
  while(!conv) {
    fak <- f(ak,...)-cc
    fbk <- f(bk,...)-cc
    ck  <- ak-(ak-bk)/(fak-fbk)*fak
    if (is.nan(ck)) {
      cat("solution non assigned \n"); return(list(solution=NULL,nit=nit))
    }
    fck <- f(ck,...)-cc
    if (fak*fck > 0) ak <- ck else bk <- ck
    conv <- max(abs(fck)) < tol
# cat(nit,ak,bk,ck,fak,fbk,fck,"\n")
    nit <- nit+1
    if (nit==maxit)  break
  }
  return(list(solution=ck,nit=nit))
}

