bagg<-function(x,y,arg,B=10,seed=1,method="worpl",propor=0.5,
estimator="greedy",M=2,m=5,splitfreq=1)
{
n<-length(y)
d<-dim(x)[2]
vals<-matrix(0,B,1)
bootstrap<-function(n,seed,method.propor){1}

if (estimator=="greedy")
for (i in 1:B){
   seed<-seed+i
   boot<-bootstrap(n,seed=seed,method=method,propor=propor)
   if (d==1) xsub<-matrix(x[boot],length(boot),1) else xsub<-x[boot,]
   ysub<-matrix(y[boot],length(boot),1)
   vals[i]<-greedy(xsub,ysub,arg,M,m,splitfreq=splitfreq)$val
}
else  # estimator=="linear"
for (i in 1:B){
   seed<-seed+i
   boot<-bootstrap(n,seed=seed,method=method,propor=propor)
   if (d==1) xsub<-matrix(x[boot],length(boot),1) else xsub<-x[boot,]
   ysub<-matrix(y[boot],length(boot),1)
   lin<-linear(xsub,ysub)
   vals[i]<-lin$beta0+lin$beta1%*%arg
}

val<-mean(vals)
return(val)
}










