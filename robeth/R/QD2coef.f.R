"QD2coef.f" <- function(lambda,yc,delta,muI,sigmaI,zero=1e-4) {
n   <- length(yc); P   <- ppoints(n)
qn  <- rep(0,n) 
tol <- 0.001; maxit <- 20
f.res <- .Fortran("qd2func",
lambda=as.double(lambda),yc=as.double(yc),delta=to.single(delta),n=to.integer(n),mui=as.double(muI),
sigmai=as.double(sigmaI),zero=as.double(zero),tol=as.double(tol),
maxit=to.integer(maxit),p=as.double(P),qn=as.double(qn))
ql   <- qloggamma(P,lambda)
z    <- lsfit(ql,f.res$qn)
mu   <- as.numeric(z$coef[1])
sig  <- as.numeric(z$coef[2])
list(lambda=lambda,mu=mu,sig=sig)}

