pp.plot<-function(dendat=NULL,compa="gauss",basis="gauss",mean=0,sig=1,df=1,
gnum=1000,d=1,R=3,pptype="1d",cex.lab=1,cex.axis=1,col="blue",lwd=1)
# basis is either data (dendat) or a theoretical distribution
{
if (pptype=="1d"){
   p<-dendat[order(dendat)]
   if (compa=="gauss") y<-pnorm(p,mean=mean,sd=sig)
   if (compa=="student") y<-pt((p-mean)/sig,df=df)
   if (compa=="unif") y<-punif((p-mean)/sig)
   if (compa=="exp") y<-pexp((p-mean)/sig)
   if (compa=="doubleexp") 
      y<-0.5*(1-pexp(-(p-mean)/sig))+0.5*pexp((p-mean)/sig)
   n<-length(dendat) #dim(dendat)[1]
   x<-seq(1:n)/n
   tyyppi<-"p"
   xlab<-"empirical distribution function"
   ylab<-"compared distribution function"
}
if (pptype=="v2p"){
      rp<-tailfunc(R,d,type=compa,gnum=gnum,sig=sig,nu=df)
      y<-rp$proba
      rp2<-tailfunc(R,d,type=basis,gnum=gnum,sig=sig,nu=df)
      x<-rp2$proba
      tyyppi="l"
      xlab<-"empirical"
      ylab<-"model"
}
if (pptype=="ddplot"){
}

plot(x,y,
type=tyyppi,
xlim=c(0,1),ylim=c(0,1),
xlab=xlab,ylab=ylab,cex.lab=cex.lab,cex.axis=cex.axis)
segments(0,0,1,1,col=col,lwd=lwd)
}

