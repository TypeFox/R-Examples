rtnegbin <-
function(N,mu,dispersion,l.bound){
sample=rep(NA,N)
mu=rep(mu,length=N)
l.bound=rep(l.bound,length=N)
n=1:N
ind=n
while(length(ind)!=0){ ##rejection sampling
sample[ind]=rnbinom(length(ind),mu=mu[ind],size=dispersion)
ind=which(sample<l.bound)
}
return(sample)
}

