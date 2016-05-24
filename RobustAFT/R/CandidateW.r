CandidateW <- function(ind,X,y,delta,MAXIT,TOL)
{
  n  <- length(y)
  p <- ncol(X)
  b  <- 0.5
  sigma0 <- 1
  s0 <- 1
  c0 <- -0.1351788
#
  options(warn=-1)
  if (p==1) sres <- survreg(Surv(y,delta)~1,     dist="extreme",subset=ind)
  else      sres <- survreg(Surv(y,delta)~X[,-1],dist="extreme",subset=ind)
  options(warn=0)
  if (is.na(sres$scale)) return(rep(NA,p))
  beta0 <- sres$coef
  sini <- s.eq.logW(X,y,N=1,delta,sigma0,b,beta0,
                 MAXIT,TOL,mu0=c0,s0=s0,ipsi=4,xk=1.717816,ialg=3,meth=4)$S
#
  rs    <- as.vector((y-X%*%as.matrix(beta0))/sini)
  indc <- (1:n)[delta==0]
  rc   <- rs[indc]
  uc   <- unique(rc)
  if (length(uc) > 0) {
    Exp.rc   <- ResExpW(uc)
    imatch   <- match(rc,uc,nomatch=0)
    rs[indc] <- Exp.rc[imatch]
  }
  cu    <- tutl(0.85)
  input <- list(tau=beta0,v=sini)
  yy    <- X%*%as.matrix(beta0) + rs*sini
  z     <- TML.logWeibull(X,yy,cu=cu,initial="input",otp="fixed",cov ="no",
         input=input,iv=1,nrep=0,seed=1313,maxit=100,tol=0.0001,gam=0.5,nitmon=FALSE)
  coef <- z$th1
  coef
}