stab.mean <-
function(Y,class,cls2=NULL,resample,times=NULL,alpha=NULL,...){
   Y=as.matrix(Y)
   tn=ncol(Y)
   if(is.null(alpha))alpha=0.05
   if(is.null(times))times=1000
   if(is.null(cls2))cls2=NULL
   if(tn==1)return(Rank2(Y[,1],class,cls2,resample,times,alpha))
   else{
      RES=list()
      for(i in 1:tn)RES[[i]]=Rank2(Y[,i],class,cls2,resample,times,alpha)
      names(RES)=colnames(Y)
      return(RES)
   }
}
