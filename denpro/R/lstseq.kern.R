lstseq.kern<-function(dendat,hseq,N,lstree=NULL,level=NULL,
Q=NULL,kernel="gauss",hw=NULL,algo="leafsfirst",support=NULL)
{
hnum<-length(hseq)
if ((hnum>1) && (hseq[1]<hseq[2])) hseq<-hseq[seq(hnum,1)]

if (algo=="leafsfirst"){

  for (i in 1:hnum){   
      h<-hseq[i]
      pcf<-pcf.kern(dendat,h,N,kernel=kernel,support=support)
      if (!is.null(lstree)) lf<-leafsfirst(pcf)
      if (!is.null(level)){ 
           lev<-level*max(pcf$value)  
           refe<-locofmax(pcf)
           st<-leafsfirst(pcf,lev=lev,refe=refe)
      }
      if (i==1){
           if (hnum==1){ 
               pcfseq<-pcf
               if (!is.null(lstree)) lstseq<-lf
               if (!is.null(level)) stseq<-st
           }
           else{
               pcfseq<-list(pcf)
               if (!is.null(lstree)) lstseq<-list(lf)
               if (!is.null(level)) stseq<-list(st)
           }
      }
      else{
          pcfseq<-c(pcfseq,list(pcf))
          if (!is.null(lstree)) lstseq<-c(lstseq,list(lf))
          if (!is.null(level)) stseq<-c(stseq,list(st))
      }
  }

}
else{  #algo=="decomdyna"
  lstseq<-profkern(dendat,hseq,N,Q,kernel=kernel,hw=hw)
}

if (is.null(lstree)) lstseq<-NULL
if (is.null(level)) stseq<-NULL
return(list(lstseq=lstseq,pcfseq=pcfseq,stseq=stseq,hseq=hseq,type="kernel"))
}

