qq.plot<-function(dendat=NULL,compa="gauss",basis="gauss",
mean=0,sig=1,df=1,
gnum=1000,d=1,R=3,qqtype="1d",cex.lab=1,cex.axis=1,col="blue",lwd=1,flip=FALSE,
xlab="compared quantiles",ylab="empirical quantiles")
{
if (qqtype=="1d"){
   n<-length(dendat) #dim(dendat)[1]
   p<-(seq(1:n)-1/2)/n
   if (compa=="gauss") x<-qnorm(p,mean=mean,sd=sig)
   if (compa=="student") x<-sig*qt(p,df=df)+mean
   if (compa=="unif") x<-sig*qunif(p)+mean
   if (compa=="exp") x<-sig*qexp(p)+mean
   if (compa=="doubleexp"){
       x<-sig*qexp(p)+mean
       alku<-which(p<0.5)
       loppu<-which(p>=0.5)
       x[alku]<--sig*qexp(1-2*p[alku])+mean
       x[loppu]<-sig*qexp(2*p[loppu]-1)+mean
   }
   y<-dendat[order(dendat)]
   tyyppi<-"p"
}
if (qqtype=="lower"){
   n<-length(dendat) #dim(dendat)[1]
   p<-(seq(1:n)-1/2)/n
   if (compa=="gauss") x<-qnorm(p/2,mean=mean,sd=sig)
   if (compa=="student") x<-sig*qt(p/2,df=df)+mean
   if (compa=="unif") x<-sig*qunif(p/2)+mean
   if (compa=="exp") x<-sig*qexp(p/2)+mean
   y<-dendat[order(dendat)]
   tyyppi<-"p"
}
if (qqtype=="p2v"){
     rp<-tailfunc(R,d,type=compa,gnum=gnum,sig=sig,nu=df)
     x<-rp$volu
     rp2<-tailfunc(R,d,type=basis,gnum=gnum,sig=sig,nu=df)
     y<-rp2$volu
     tyyppi="l"
     ylab<-"empirical"
     xlab<-"model"
}

if (!flip){
plot(x,y,type=tyyppi,ylab=ylab,xlab=xlab,cex.lab=cex.lab,cex.axis=cex.axis)
maxxy<-max(max(x),max(y))
minxy<-min(min(x),min(y))
segments(minxy,minxy,maxxy,maxxy,col=col,lwd=lwd)
}
if (flip){
 plot(y,x,type=tyyppi,ylab=xlab,xlab=ylab,cex.lab=cex.lab,cex.axis=cex.axis)
 maxxy<-max(max(x),max(y))
 minxy<-min(min(x),min(y))
 segments(minxy,minxy,maxxy,maxxy,col=col,lwd=lwd)
}

}





