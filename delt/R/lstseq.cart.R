lstseq.cart<-function(treeseq,maxleaf=NULL,lstree=NULL,level=NULL,indvec=NULL)
{

if ((is.null(indvec)) && (is.null(maxleaf))){ 
    leaf<-treeseq$leafs
    alfa<-treeseq$alfa
    alkm<-length(leaf)
}
else if (!is.null(maxleaf)){
    aplnum<-roundlnum(treeseq$leafs,maxleaf)
    indeksi<-detsi(treeseq$leafs,aplnum)
    leaf<-treeseq$leafs[indeksi:length(treeseq$leafs)]
    alfa<-treeseq$alfa[indeksi:length(treeseq$leafs)]
    alkm<-length(leaf)
}
else if (!is.null(indvec)){ 
  leaf<-treeseq$leafs[indvec]
  alfa<-treeseq$alfa[indvec]
  alkm<-length(indvec)
}

tuloleaf<-matrix(0,alkm,1)
tuloalfa<-matrix(0,alkm,1)
laskuri<-1

for (inds in alkm:1){  # start with the oversmoothed estimate 
     leafnum<-leaf[inds]
   
     pcf<-eval.pick(treeseq,leafnum)
     #pv<-partition(pcf,suppo)
     #lst<-profgene(pv$values,pv$recs,frekv=F,cvol=T,ccen=T,cfre=F)
     #lst<-proftree(pcf)
     if (!is.null(lstree)) lst<-leafsfirst(pcf)
     if (!is.null(level)){ 
           lev<-level*max(pcf$value)  
           refe<-locofmax(pcf)
           st<-leafsfirst(pcf,lev=lev,refe=refe)
     }

     if (inds==alkm){
        if (alkm==1){
               pcfseq<-pcf
               if (!is.null(lstree)) lstseq<-lst
               if (!is.null(level)) stseq<-st
       }
        else{
               pcfseq<-list(pcf)
               if (!is.null(lstree)) lstseq<-list(lst)
               if (!is.null(level)) stseq<-list(st)  
        }
     }
     else{
          pcfseq<-c(pcfseq,list(pcf))
          if (!is.null(lstree)) lstseq<-c(lstseq,list(lst))
          if (!is.null(level)) stseq<-c(stseq,list(st))
     }
     
    tuloleaf[laskuri]<-leaf[inds]
    tuloalfa[laskuri]<-alfa[inds]
    laskuri<-laskuri+1
}

if (is.null(lstree))  lstseq<-NULL
if (is.null(level)) stseq<-NULL
return(list(lstseq=lstseq,pcfseq=pcfseq,stseq=stseq,hseq=tuloalfa,
leaf=tuloleaf,type="carthisto"))
}

