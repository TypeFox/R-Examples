pgram <-function(x, fftlen,...){

 pgram_full<-function(x,fftlen,np1,np2,halflen,alpha,rejalpha,logsw,datastr,typeci,typepgram, colci,colpgram)
 {
    nx=length(x)

    if (np1 <= halflen)

    {stop("'np1' is too small")}

    x[is.nan(x)]=0 
    nav=floor(nx/fftlen)                         
    start=1
    stop=start+fftlen-1
    xtsum<-matrix(1,fftlen,1)


    for (j in 1:nav)  { 
        xin=x[start:stop]
xt=abs(fft(xin))        
xtsum=xtsum+(abs(xt)^2)/(2*pi*fftlen)    
 start=start+fftlen
stop=stop+fftlen
      }

     xtav=xtsum/nav       

     CH1=qchisq(1-rejalpha,2*nav)/qchisq(.5,2*nav)  
                                              
     frat<-matrix(0,1,np2-np1+1)
     backav<-matrix(0,1,np2-np1+1)
     fthr<-matrix(0,1,np2-np1+1)

    for (j in np1:np2) {
       index=c(seq((j-halflen),(j-1)),seq((j+1),(j+halflen)))
        medback=median(xtsum[index])/(2*nav)             
        rejthr=medback*CH1
        good<-which( (xtsum/(2*nav)) < rejthr )
        numgood=length(good)
       
        back= sum(xtsum[good])/(2*nav*numgood)   
        num=xtsum[j]/(2*nav) 
        frat[j]=num/back
        backav[j]=back
        fthr[j]=qf(1-alpha,2*nav,2*nav*numgood)
      }  

    if (logsw)  {

       indexp=c(np1:np2)
       pgrm=fthr[indexp]*backav[indexp]*(2*nav)/nav
       matplot(indexp, log10(xtav[indexp]),xlab="frequency index j", ylab="periodogram (log10 scale)", type="l", lwd=1, lty=4, col="blue")
       lines(indexp, log10(pgrm), type="l", lwd=1, lty=4, col="red")
       title(main=paste("Periodogram based on FFT of length:",fftlen), sub=(paste('the nr of FFT: ',nav,'; alpha: ',alpha,'; half of smoothing  window: m = ',halflen) ))
      } else {
       indexp=c(np1:np2)
       pgrm=fthr[indexp]*backav[indexp]*(2*nav)/nav
       plot(indexp, xtav[indexp],xlab="frequency index j", ylab="periodogram ", type=typeci, lwd=1, col=colci)
       lines(indexp, pgrm,  type=typepgram, lwd=1, col=colpgram)
       title(main=paste("Periodogram based on FFT of length:",fftlen), sub=(paste('the nr of FFT: ',nav,'; alpha: ',alpha,'; half of smoothing  window: m = ',halflen) ))
     }
                                                
      fmax=max(frat)
      imax<-match(max(frat), frat)   
      xtpv=1-pf(fmax,2*nav,2*nav*numgood)
 
      cat(paste('Most significant periodogram bin ',imax-1,' P-value : ',xtpv,'\n'))     

}
 
L<-modifyList(list(np1=5,np2=fftlen/2,halflen=4, alpha=0.05, rejalpha=0.01,logsw=1,datastr="data",typeci="b",typepgram="b", colci="red",colpgram="red"), list(x = x, fftlen=fftlen,...))

 do.call(pgram_full,L)

 }


