rosunc <-
function (modeltype,w,s,delta,mx.dead,h,m,u,slope,sdu=0,sdm=0,sds=0,sdw=0,sdd=0,nsim=1000) {

# correlations structure of moisture values
#ms<-data(scenarios)
#ms5<-data.frame(ms[,c(8,9,10,4,3)])
#round(cor(ms5),3)
#sig = matrix(c(1,1,.99,.1,.1,1,1,.99,.1,.1,.99,.99,1,.08,.08,.1,.1,.08,1,1,.1,.1,.08,1,1),nrow=5)
#library(mvtnorm)
#m<-rmvnorm(1000,mean=as.numeric(ex[,11:15]),sigma=sig,method="svd") 

# sampling from gaussian distribution 
u<-rnorm(nsim,u,u*sdu)
if (is.data.frame(m)==F) m<-t(data.frame(m))
if (is.data.frame(w)==F) w<-t(data.frame(w))

m0<-which(m[,1:5]==0)
w0<-which(w[,1:5]==0)

m<-data.frame(
  rnorm(nsim,m[,1],m[,1]*sdm),
  rnorm(nsim,m[,2],m[,2]*sdm),
  rnorm(nsim,m[,3],m[,3]*sdm),
  rnorm(nsim,m[,4],m[,4]*sdm),
  rnorm(nsim,m[,5],m[,5]*sdm))
m[,m0]<-0

if (length(which(m[,4]>0 & m[,4]<30))>0) 
  {m[which(m[,4]>0 & m[,4]<30),4]<-30}
slope<-rnorm(nsim,slope,slope*sds)
w<-data.frame(
  rnorm(nsim,w[,1],w[,1]*sdw),
  rnorm(nsim,w[,2],w[,2]*sdw),
  rnorm(nsim,w[,3],w[,3]*sdw),
  rnorm(nsim,w[,4],w[,4]*sdw),
  rnorm(nsim,w[,5],w[,5]*sdw))
w[,w0]<-0
delta<-rnorm(nsim,delta,delta*sdd)

pred<-numeric()
for (n in 1:nsim) {
  pred[n]<-ros(modeltype,w[n,],s,delta[n],mx.dead,h,m[n,],u[n],slope[n])[15]
}

return (as.numeric(pred))

}
