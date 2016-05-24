cm.n.clopper.pearson <-
function(p,size,cm.effect,alpha=0.1,uniroot.lower=k+1,uniroot.upper=1e+100,uniroot.tol=1e-10,uniroot.maxiter=100000){
k=sum(size)
if(alpha<0 || alpha>1) {
stop("'alpha' must be a number between 0 and 1")}
if(p<0 || p>1) {
stop("'p' must be a number between 0 and 1")}
K<-c(0:k)
xi<-dgbinom(K,size,cm.effect)
xi<-rev(xi)
nstar<-uniroot(function(n) xi%*%pbeta(p,K+1,n-K)-(1-alpha),lower=uniroot.lower,upper=uniroot.upper,tol=uniroot.tol,maxiter=uniroot.maxiter)
return(ceiling(nstar$root))
}
