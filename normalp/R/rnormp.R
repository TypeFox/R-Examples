rnormp<-function(n,mu=0,sigmap=1,p=2,method=c("def","chiodi")){
if(!is.numeric(n)||!is.numeric(mu)||!is.numeric(sigmap)||!is.numeric(p)) 
stop (" Non-numeric argument to mathematical function")
if(p<1) stop("p must be at least equal to one")
if(sigmap<=0) stop("sigmap must be positive")
method <- match.arg(method)
if (method=="def"){
qg<-rgamma(n,shape=1/p,scale=p)
z<-qg^(1/p)
z<-ifelse(runif(n)<0.5,-z,z)
x<-mu+z*sigmap
}
if (method=="chiodi"){
i<-0
x<-c(rep(0,n))
while (i<n){
u<-runif(1,-1,1)
v<-runif(1,-1,1)
z<-abs(u)^p+abs(v)^(p/(p-1))
if (z<=1){
i<-i+1
x[i]<-u*(-p*log(z)/z)^(1/p)
x[i]<-mu+x[i]*sigmap
}
}
}
x
}

