slbind.cond<-function(intvar, statslist, var.suffix, sl.ind=NULL, who.evs=NULL, type=1, ...){
   if(type==1&&length(intvar)!=length(statslist)){
   stop("The interaction variable must be the same length as the statslist for global statistics.")
   }

   if(!is.null(who.evs)&&type==1){
   stop("The who.evs parameter is only valid for local statistics.")
   }

   if(is.null(sl.ind)){
   warning("No sufficient statistic column indices were supplied. Applying the interaction term to every (non-reference) variable.")
   sl.ind<-1:ncol(statslist[[1]][[type]][1,,])
   }

   if(is.null(who.evs)){
   who.evs<-1:length(statslist)
   }


   sl.copy<-statslist

   for(i in who.evs){
      for(j in sl.ind){
         bplt<-statslist[[i]][[type]][,,j]
         bplt[bplt==1]<-intvar[i]
         sl.copy[[i]][[type]][,,j]<-bplt
      }
       bnames<-dimnames(sl.copy[[i]][[type]])
       bnames[[3]]<-c(bnames[[3]],paste(dimnames(statslist[[i]][[type]])[[3]][sl.ind],var.suffix,sep="."))
      sl.copy[[i]][[1]]<-abind(statslist[[i]][[1]],sl.copy[[i]][[1]][,,sl.ind],new.names=bnames,...)
   }
   sl.copy
}
