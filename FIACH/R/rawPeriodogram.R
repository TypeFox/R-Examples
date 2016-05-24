rawPeriodogram<-function(x,sf=NULL,plot=TRUE,amp=FALSE,phase=FALSE,N=NULL){
  x<-as.matrix(x)
  n<-nrow(x)
  
  if(is.null(N)){N<-nextn(n)}
  if(ncol(x)>1 & plot==TRUE){warning("Number of series is greater than 1: Only the first is plotted")}
  
  nyq.ind<-floor(N/2+1)
  nyq.ind
  fft.sig<-fftN(x,N)
  ampl<-Mod(fft.sig)
  if(phase==TRUE){pha<-Arg(fft.sig)}
  pow<-(ampl^2)/floor(N/2)
  pow.per<-pow[2:nyq.ind,,drop=FALSE]

  
  if(is.null(sf)){freq<-(1:(floor(N/2)))/(N)}else{
    freq<-((1:(floor(N/2)))/(N))*sf
  }
  if(is.null(sf)){xlab<-"Normalised Frequency (Hz)"}else{
    xlab<-"Frequency (Hz)"
  }
  
  if(plot==TRUE){plot(freq,pow.per[,1],type="h",lwd=2,ylab="Power",xlab=xlab)}
  ret<-list()
  ret$power<-pow.per
  if(amp){ret$amplitude<-ampl[2:nyq.ind,,drop=FALSE]}
  if(phase){ret$phase<-pha[2:nyq.ind,,drop=FALSE]}
  ret$freq<-freq
  class(ret)<-"spectrum"
  return(ret)
}
