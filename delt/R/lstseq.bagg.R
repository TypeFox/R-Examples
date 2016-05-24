lstseq.bagg<-function(dendat,B,
lstree=NULL,level=NULL,
maxleaf=NULL,leafseq=NULL,
minobs=NULL,seed=1,sample="bagg",prune="off",
splitscan=0,seedf=1,scatter=0,src="c",method="loglik")
{
n<-dim(dendat)[1]
if (is.null(minobs)) minobs<-ceiling(sqrt(n)/2)

if (!is.null(maxleaf)){  
  leafseq<-seq(maxleaf,2)
  hnum<-maxleaf-1
}
else{
  hnum<-length(leafseq)
  if ((hnum>1) && (leafseq[1]>leafseq[2])) leafseq<-leafseq[seq(hnum,1)]
}

for (i in 1:hnum){   
      leaf<-leafseq[i]
      pcf<-eval.bagg(dendat,B,leaf,
           minobs=minobs,seed=seed,
           sample=sample,prune=prune,
           splitscan=splitscan,seedf=seedf,
           scatter=scatter,src=src,method=method)


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
return(list(lstseq=lstseq,pcfseq=pcfseq,stseq=stseq,
hseq=leafseq,type="bagghisto"))
}

