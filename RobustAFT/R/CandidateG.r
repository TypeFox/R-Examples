CandidateG <-
function(ind,X,y,delta,MAXIT,TOL) {
n  <- length(y); p <- ncol(X); b  <- 0.5; sigma0 <- 1; s0 <- 1
#
options(warn=-1)
beta0 <- lsfit(X[ind,],y[ind],intercept=FALSE)$coef
options(warn=0)
sini  <- s.eq.Gauss(X,y,N=1,delta,sigma0,b,beta0,
         MAXIT,TOL,s0=s0,ipsi=4,xk=1.5477,lint=1,ialg=3,meth=4)$S
#
rs   <- as.vector((y-X%*%as.matrix(beta0))/sini)
indc <- (1:n)[delta==0]
rc   <- rs[indc]
uc   <- unique(rc)
if (length(uc) > 0) {
  Exp.rc   <- ResExpG(uc)
  imatch   <- match(rc,uc,nomatch=0)
  rs[indc] <- Exp.rc[imatch] }
cu    <- 0.6745
input <- list(lambda=beta0,sigma=sini)
yy    <- X%*%as.matrix(beta0) + rs*sini
z     <- TML.gauss(X,yy,cu=cu,initial="input",otp="fixed",cov ="no",input=input,
         iv=1,nrep=0,tol=0.0001,seed=1313)
coef  <- z$th1
coef}

