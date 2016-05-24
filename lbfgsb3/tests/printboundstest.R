rm(list=ls())
require("lbfgsb3")
#####################
# Simple bounds test (use n=4)

bt.f<-function(x){
 sum(x*x)
}

bt.g<-function(x){
  gg<-2.0*x
}

n<-4
lower<-rep(0,n)
upper<-lower # to get arrays set
# bdmsk<-rep(1,n)
# bdmsk[(trunc(n/2)+1)]<-0
for (i in 1:n) { 
    lower[i]<-1.0*(i-1)*(n-1)/n
    upper[i]<-1.0*i*(n+1)/n
}
xx<-0.5*(lower+upper)

abt<-lbfgsb3(xx, bt.f, bt.g, lower=lower, upper=upper, control=list(trace=1, iprint=200L))
print(abt)

