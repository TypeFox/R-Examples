test3=function(Nsim=10^4,lambda=100){
  spread=3*sqrt(lambda)
  t=round(seq(max(0,lambda-spread),lambda+spread,1))
  prob=ppois(t,lambda)
  X=rep(0,Nsim)
  for (i in 1:Nsim){
     u=runif(1)
     X[i]=t[1]+sum(prob<u)-1 }
  }

test4=function(Nsim=10^4,lambda=100){
 rpois(Nsim,lambda)}

