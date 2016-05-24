rf2tree<-function(rf,support,N=NULL)
{
forest<-rf$forest
d<-length(support)/2

if (is.null(N)) N<-rep(60,d)

nr<-length(forest$ndbigtree)     #number of trees in the forest
glob<-1
while (glob <= nr){

   # build trnew

   left<-forest$leftDaughter[,glob]
   right<-forest$rightDaughter[,glob]
   direc<-forest$bestvar[,glob]
   direc[(forest$bestvar[,glob]==0)]<-NA
   mean<-forest$nodepred[,glob]
   nodenum<-length(forest$xbestsplit[,glob])
   nelem<-rep(1,nodenum)
   split<-matrix(NA,nodenum,1)
   for (i in 1:nodenum){
       vec<-direc[i]
       ala<-support[2*vec-1]
       yla<-support[2*vec]
       splitti<-round(N[vec]*(forest$xbestsplit[i,glob]-ala)/(yla-ala))
       split[i]<-min(max(splitti,1),(N[vec]-1))
   }
   split[(forest$bestvar[,glob]==0)]<-NA

   trnew<-list(split=split,direc=direc,mean=mean,nelem=nelem,
   left=left,right=right,#low=low,upp=upp,
   N=N,support=support)

   lu<-lowupp(trnew)
   trnew$low<-lu$low
   trnew$upp<-lu$upp
   trnew$volume<-rep(1,nodenum)

   if (glob>1) tr<-treeadd(tr,trnew) else tr<-trnew

   glob<-glob+1
}

dh<-downhigh(tr)
tr$down<-dh$down
tr$high<-dh$high   
tr$value<-dh$value
tr$infopointer<-dh$infopointer

return(tr)
}
