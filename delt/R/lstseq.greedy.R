lstseq.greedy<-function(dendat,maxleaf,lstree=NULL,level=NULL)
#treeseq,lsets=FALSE,invalue=FALSE,indvec=NULL,
{
hseq<-seq(maxleaf,1)
hnum<-length(hseq)

for (i in 1:hnum){

     leaf<-hseq[i]
     pcf<-eval.greedy(dendat,leaf=leaf)

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

if (is.null(lstree))  lstseq<-NULL
if (is.null(level)) stseq<-NULL
return(list(lstseq=lstseq,pcfseq=pcfseq,stseq=stseq,hseq=hseq,type="greedy"))
}




