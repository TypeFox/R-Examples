n.clopper.pearson <-
function(k,p,alpha=0.1,uniroot.lower=k+1,uniroot.upper=1e+100,uniroot.maxiter=100000,uniroot.tol=1e-10){
l<-round(k)
if (is.na(k) || k < 0 || max(abs(k - l)) > 1e-07) 
        stop("'k' must be nonnegative and integer")
if (p<=0 || p>=1 ){
stop("'p'  must be a number between 0 and 1")}
if(alpha<0||alpha>=1){
stop("'alpha'  must be a number between 0 and 1")}
nstar<-uniroot(function(n) pbeta(p,k+1,n-k)-(1-alpha),lower=uniroot.lower,upper=uniroot.upper,tol=uniroot.tol,maxiter=uniroot.maxiter)
return(ceiling(nstar$root))}
