
N<-c(50,100,200)
nsim<-length(N)
stout<-matrix(NA,nsim,3,dimnames=list(paste("n=",N,sep=""),c("simple","EMICM","computeMLE")))
set.seed(1)


for (i in 1:nsim){
n<-N[i]
l<-exp(rnorm(n))
r<-l+exp(rnorm(n))
g<-c(rep(0,n/2),rep(1,n/2))
st0<-system.time(out<-ictest(l,r,g,initfit=NULL,icontrol=icfitControl(epsilon=1e-6)))
st1<-system.time(out<-ictest(l,r,g,initfit="initEMICM",icontrol=icfitControl(epsilon=1e-6)))
st2<-system.time(out<-ictest(l,r,g,initfit="initcomputeMLE",icontrol=icfitControl(epsilon=1e-6)))
stout[i,1]<-st0[1]
stout[i,2]<-st1[1]
stout[i,3]<-st2[1]
}

stout