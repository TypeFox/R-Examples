kmms <-
function(dd,gam=0.95){
#    Kaplan- Meier(KM) estimate of mean and Standard Error of
#    the mean for left censored data  see Section 3.5

t1<- plend(dd)

t1<- t1[rev(order(t1[,1])),]
x<- t1[,1] ; n<-t1$n  ; d<-t1$r
ple<- c(t1$ple[2:(dim(t1)[1])],0)
sv<- 1 - ple
xv<- abs( diff(c(t1[,1],0) ) )
a<- sv*xv

ac<-xv*ple
B<- rev(cumsum(rev(ac)))

t2<-data.frame(x,n,d,ple,sv,a,ac,B,xv)
t2<- t2[order(t2$x),]
nx <- sum(dd[,2])  #  CORRECTION 14 JUN 2006 nx<- dim(t2)[1]
#  compute terms for variance
vb<-ifelse( (t2[,2] -t2[,3])==0 , 0, 1/((t2[,2] -t2[,3])*t2[,2]) )
vb<- t2$d*t2$B^2*vb
t2<-data.frame(t2,vb)  #  t.2<<-t2
kmm<-sum( t2$a) # Kaplan-Meier mean
kmvb<-sum(vb)   # Kaplan-Meier variance (unadjusted)
kmseu<- sqrt(kmvb); kmse<- kmseu* sqrt(nx/(nx-1) )
# 
cd<- qt( gam ,nx-1)*kmse
kml<- kmm - cd ; kmu<- kmm + cd
out<-list("KM.mean"=kmm,"KM.LCL"=kml,"KM.UCL"=kmu,"KM.se"=kmse,"gamma"=gam)
out
}

