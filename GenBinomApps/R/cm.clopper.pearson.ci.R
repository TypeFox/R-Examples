cm.clopper.pearson.ci <-
function(n,size,cm.effect,alpha=0.1,
CI="upper",uniroot.lower=0,uniroot.upper=1,uniroot.maxiter=100000,uniroot.tol=1e-10){
k=sum(size)
l<-round(size)
if (any(is.na(size) | (size < 0)) || max(abs(size - l)) > 1e-07) 
        stop("'size' must be nonnegative and integer")
m<-round(n)
if (is.na(n) || n < k || max(abs(n - m)) > 1e-07) 
        stop("'n' must be nonnegative and integer >= k")
if(alpha<0 || alpha>1) {
stop("'alpha' must be a number between 0 and 1")}

K<-c(0:k)
xi<-dgbinom(K,size,cm.effect)
xi<-rev(xi)

if (CI=="upper") 
{ll<-0
 ul<-uniroot(function(pii) xi%*%pbeta(pii,K+1,n-K)-(1-alpha),lower=uniroot.lower,upper=uniroot.upper,tol=uniroot.tol,maxiter=uniroot.maxiter)$root
}
else if (CI=="lower")
{ul<-1
if (xi[1]>=(1-alpha)){
ll<-0}
else{
K <- c(1:k)
        ll <- uniroot(function(pii) xi %*% c(pbeta(pii,1e-100,1+n), pbeta(pii, K, 1 + 
            n - K)) - alpha, lower = uniroot.lower, upper = uniroot.upper, 
            tol = uniroot.tol, maxiter = uniroot.maxiter)$root
}
}
else if (CI=="two.sided")
{
if (xi[1]>=(1-alpha/2)){
ll<-0}
else{
 Kl <- c(1:k)
        ll <- uniroot(function(pii) xi %*% c(pbeta(pii,1e-100,1+n), pbeta(pii, Kl, 
            1 + n - Kl)) - alpha/2, lower = uniroot.lower, upper = uniroot.upper, 
            tol = uniroot.tol, maxiter = uniroot.maxiter)$root
}
ul<-uniroot(function(pii) xi%*%pbeta(pii,K+1,n-K)-(1-alpha/2),lower=uniroot.lower,upper=uniroot.upper,tol=uniroot.tol,maxiter=uniroot.maxiter)$root
}
else stop("undefined CI detected")
data.frame(Confidence.Interval=CI,Lower.limit=ll,Upper.limit=ul,alpha=alpha,row.names="")
}
