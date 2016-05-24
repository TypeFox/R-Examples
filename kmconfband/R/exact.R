exact<-function(sobj,x){

n<-length(sobj$time)

select<-(1:n)*as.numeric(sobj$n.event>0)

k<-sum(sobj$n.event>0)

bounds<-matrix(0,k+1,2)

K<-function(phat,y) phat*log(phat/y)+(1-phat)*log((1-phat)/(1-y))

cvn<-x

bounds[1,]<-c(exp(-cvn),1)

for (j in seq(along=1:(k-1))){

usubi<-(sobj$surv[select])[j]


f<-function(p) K(usubi,p)-cvn

xin<-zbrent(f,c(1.0e-80,usubi),1.0e-10)

bounds[j+1,1]<-xin

xip<-zbrent(f,c(usubi,1.0-1.0e-80),1.0e-10)

bounds[j+1,2]<-xip}

usubi<-(sobj$surv[select])[k]


f<-function(p) K(usubi,p)-cvn

if (usubi==0) bounds[k+1,]<-c(0,1-exp(-cvn))
	else bounds[k+1,]<-c(zbrent(f,c(0.00000000001,usubi),1.0e-10),zbrent(f,c(usubi,0.99999999999),1.0e-10))


bounds}
