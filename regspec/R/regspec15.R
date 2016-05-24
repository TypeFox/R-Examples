regspec<-function(D,deltat=1,nb=100,varmult=2,smthpar=0.8,ebeta=NULL,
	vbeta=NULL,filter=NULL, freq.out=seq(0,0.5,length=200),
	plot.spec=TRUE, plot.log=FALSE, plot.pgram=FALSE,
	plot.intervals=TRUE, ylim=NULL, 
	SARIMA=list(sigma2=1),centred=FALSE,intname=NULL,
	xlab="Frequency", ylab="Spectrum", ...){

#
# Check length of data in "starting" case
#
nD<-length(D)

if (is.null(ebeta) && is.null(vbeta) & nD < 2)
	stop("With no prior data, length of series has be greater than 2")

if(is.null(ebeta)){ebeta<-rep(0,nb)}

if(is.null(vbeta)){
vbeta<-vbetafun(smthpar,nb,varmult)
}



nw<-ceiling((nD+1)/2)
fftx<-fft(D)
spec<-Mod(fftx[1:nw])^2/nD
w<-(seq(-1,1,length=nw)*0.99/4+0.25)/deltat
logI<-log(spec*deltat)
nyq<-(2*deltat)^-1

nb<-length(ebeta)
b0<-ebeta

w1<-outer(w,(1:ceiling(deltat/2)-1)*2*nyq,FUN="+")
if(deltat>1){w2<-outer(-w,(1:floor(deltat/2))*2*nyq,FUN="+")}else{w2<-NULL}
wall<-cbind(w1,w2)

#################
################# compute filter gain
sgain<-matrix(1,nrow(wall),ncol(wall))
if(!is.null(filter)){
for(i in 1:ncol(wall)){
args<--2*pi*outer(1:length(filter)-1,wall[,i])
expons<-complex(length(args),modulus=1,argument=args)
expons<-matrix(expons,length(filter),)
sgain[,i]<-(Mod(filter%*%expons))^2}
}

##################
##################

# iterate linearization
b0<-ebeta
pbeta<-solve(vbeta)
alphas<-rep(1,nw)


#for(linit in 1:30){

f0<-wall
for(i in 1:ncol(wall)){f0[,i]<-exp(basis(wall[,i],nb)%*%b0)}
f0<-f0*matrix(SARIMAspec(SARIMA,freq=c(wall))$spec,nrow(f0),ncol(f0))*sgain

s0<-rowSums(f0)
B<-matrix(0,nw,nb)
for(i in 1:ncol(wall)){B<-B+sweep(basis(wall[,i],nb),1,f0[,i]/s0,FUN="*")}

eerr<-rep(digamma(1),nw)
verr<-rep(trigamma(1),nw)
if(nD%%2==0){eerr[c(1,nw)]<-digamma(0.5)+log(2);verr[c(1,nw)]<-trigamma(0.5)}
if(!nD%%2==0){eerr[1]<-digamma(0.5)+log(2);verr[1]<-trigamma(0.5)}
if(centred==TRUE){verr[1]<-999;logI[1]<-0}

V<-B%*%vbeta%*%t(B)+diag(verr)
Vinv<-solve(V)
C<-vbeta%*%t(B)
ebeta<-ebeta+C%*%Vinv%*%(logI-log(s0)-eerr)
vbeta<-vbeta-C%*%Vinv%*%t(C)

#}

B<-basis(freq.out,nb)
g0<-log(SARIMAspec(SARIMA,freq=freq.out)$spec)
gmean<-B%*%ebeta+g0
gvar<-colSums(t(B)*vbeta%*%t(B))

vareig<-eigen(vbeta)

specmean<-exp(gmean+0.5*gvar)

interval90<-matrix(,length(freq.out),2)
interval50<-matrix(,length(freq.out),2)
for(i in 1:length(freq.out)){
interval90[i,]<-Gaussbound(gmean[i],gvar[i]^0.5,.95)
interval50[i,]<-Gaussbound(gmean[i],gvar[i]^0.5,.5)
}
colnames(interval90)<-c("lower","upper")
colnames(interval50)<-c("lower","upper")

# This is a mess!
if(plot.spec==TRUE){

  if(plot.log==TRUE){
    if(is.null(ylim)){ylim<-range(interval90)}

    par(mfcol=c(1,1),oma=c(3,0,0,0))
    plot(NA,NA,type="l",col=rgb(0.2,0.5,0.7),lwd=2,ylim=ylim,xlab=xlab,
	ylab=ylab, xlim=c(0,0.5),xaxp=c(0,0.5,5), ...)
    mtext("frequency",side=1,line=2,cex=par()$cex)
    axfreq<-round(seq(0,0.5,length=6),1)
    axis(1,line=4,at=axfreq,labels=paste(round(axfreq^-1,1)))
    mtext(paste("wavelength",intname,sep=","),side=1,line=6,cex=par()$cex)

    if(plot.intervals==TRUE){
      polygon(c(freq.out,rev(freq.out)),c(interval90[,1],rev(interval90[,2])),col=rgb(.2^.2,.5^.2,.7^.2),border=NA)
      polygon(c(freq.out,rev(freq.out)),c(interval50[,1],rev(interval50[,2])),col=rgb(.2^.5,.5^.5,.7^.5),border=NA)
      }

    points(freq.out,gmean,type="l",col=rgb(0.2,0.5,0.7),lwd=2)

    if(plot.pgram==TRUE){points(wall,log(rep(spec,deltat)),lwd=2,col=rgb(1-nw^-0.1,1-nw^-0.1,1-nw^-0.1))}
  }

if(plot.log==FALSE){
  if(is.null(ylim)){ylim<-range(exp(interval90))}
  par(mfcol=c(1,1),oma=c(3,0,0,0))
  plot(NA,NA,type="l",col=rgb(0.2,0.5,0.7),lwd=2,ylim=ylim,xlim=c(0,0.5),
	xlab=xlab,ylab=ylab,xaxp=c(0,0.5,5), ...)
  mtext("frequency",side=1,line=1,cex=par()$cex)
  axfreq<-round(seq(0,0.5,length=6),1)
  axis(1,line=4,at=axfreq,labels=paste(round(axfreq^-1,1)))
  #mtext("wavelength",side=1,line=5,cex=par()$cex)
  mtext(paste("wavelength",intname,sep=","),side=1,line=5,cex=par()$cex)

  if(plot.intervals==TRUE){
    polygon(c(freq.out,rev(freq.out)),exp(c(interval90[,1],rev(interval90[,2]))),col=rgb(.2^.2,.5^.2,.7^.2),border=NA)
    polygon(c(freq.out,rev(freq.out)),exp(c(interval50[,1],rev(interval50[,2]))),col=rgb(.2^.5,.5^.5,.7^.5),border=NA)
  }

  points(freq.out,exp(gmean),type="l",col=rgb(0.2,0.5,0.7),lwd=2)

  if(plot.pgram==TRUE){points(wall,rep(spec,deltat),lwd=2,col=rgb(1-nw^-0.1,1-nw^-0.1,1-nw^-0.1))}
}

}
#mfcol=c(1,1)
return(list(freq=freq.out,spec=specmean,logspec=gmean,interval=exp(interval90),ebeta=ebeta,vbeta=vbeta,pgram=list(freq=w,spec=spec)))
}
