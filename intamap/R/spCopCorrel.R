`covar` <-
function(h,delta,corrmodel){

covariance<-switch(corrmodel,
Gau=(1-delta[1])*(exp((-h^2)/(delta[2]^2))),
Exp=(1-delta[1])*(exp((-h)/(delta[2]))),
Sph=(1-delta[1])*(1-(3*h/(2*delta[2])-h^3/(2*delta[2]^3))),
#Mat=if(delta[3]<20,(1-delta[1])*((2^(1-delta[3]))/(gamma(delta[3]))*(h^delta[3])/(delta[2]^delta[3])*besselK(h/delta[2],delta[3]))
#Mat=(1-delta[1])/(2^(delta[3] - 1)*gamma(delta[3]))*(2*delta[3]^(1/2)*h/delta[2])^delta[3]*besselK(2*delta[3]^(1/2)*h/delta[2],delta[3])
Ste=maternmodel(h,delta)
)
if(corrmodel=="Sph"){
covariance[h>delta[2]]<-0
}
covariance[row(covariance)==col(covariance)]<-1
covariance
}

#`maternmodel` <-
#function(h,delta){
#matern<-besselK(2*delta[3]^(1/2)*h/delta[2],delta[3])
#multipl<-(1-delta[1])/(2^(delta[3] - 1)*gamma(delta[3]))*(2*delta[3]^(1/2)*h/delta[2])^delta[3]
#ifelse(matern==0 | !is.finite(multipl),0,ifelse(!is.finite(matern),1-delta[1],
#multipl*matern))
#}

`maternmodel` <-
function(h,delta){
hnew=h[lower.tri(h)]
matern<-besselK(2*delta[3]^(1/2)*hnew/delta[2],delta[3])
multipl<-(1-delta[1])/(2^(delta[3] - 1)*gamma(delta[3]))*(2*delta[3]^(1/2)*hnew/delta[2])^delta[3]
covnew<-ifelse(matern==0 | !is.finite(multipl),0,ifelse(!is.finite(matern),1-delta[1],
multipl*matern))
cov=matrix(0,nrow=dim(h)[1],ncol=dim(h)[2])
cov[lower.tri(h)]<-covnew
cov+t(cov)
}

