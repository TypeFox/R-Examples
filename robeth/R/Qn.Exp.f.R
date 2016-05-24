"Qn.Exp.f" <- function(p,yc,delta,mu,sigma,lambda,zero=1e-4) {
n  <-length(yc)
a <- -10; b <- 6; tol <- 0.001; maxit <- 20
f.res <- .Fortran("qnexp",
p=as.double(p),yc=as.double(yc),delta=to.single(delta),n=to.integer(n),mu=as.double(mu),
sigma=as.double(sigma),lambda=as.double(lambda),zero=as.double(zero),a=as.double(a),
b=as.double(b),tol=as.double(tol),maxit=to.integer(maxit),qj=double(1),itr=integer(1),
iterm=integer(1))
f.res$qj
}


