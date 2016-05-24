"Fn.Exp.f" <- function(z,y,delta,mu,sigma,lambda,zero=1e-4) {
n  <-length(y)
f.res <- .Fortran("fnexp",
z=as.double(z),y=as.double(y),delta=to.single(delta),n=to.integer(n),mu=as.double(mu),
sigma=as.double(sigma),lambda=as.double(lambda),zero=as.double(zero),res=double(1))
f.res$res
}


