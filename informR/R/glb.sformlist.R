glb.sformlist<-function(evl,sforms,new.names,dichot=TRUE,cond=FALSE,interval=FALSE,warn=TRUE){
   if(length(new.names)!=length(sforms)){
   stop("The extent of sforms must be the same as new.names.")
   }
   if(any(unlist(lapply(sforms,length))<2)){
   stop("An element of sforms is too small (sform<2).")
   }

   beta.out<-lapply(sforms,gen.sformlist,evl=evl,cond=cond,interval=interval,warn=warn)
   names(beta.out)<-new.names
   #A convenience function for iterating over the 3rd element of each statslist array
   lcomb<-function(x,new.name="X"){
     tmp<-x[,,1] 
     for(i in 2:dim(x)[[3]]){
       tmp<-x[,,i]+tmp
       if(dichot && any(tmp>1) ){
       tmp[which(tmp>1)]<-1}
      }
      outpm<-array(tmp,dim=c(dim(x)[1:2],1),dimnames=list(dimnames(x)[[1]],dimnames(x)[[2]],new.name))
      return(outpm)
   }
   #a convenience function for inheriting names from abind via do.call
   gabind<-function(...) abind(...,new.names=list(NULL,NULL,new.names))

   lb1<-lapply(beta.out,lapply,lcomb)

   uids<-names(evl$eventlist)
   tmp<-list()
   for(i in 1:length(uids)){
      tmp[[i]]<-lapply(lb1,function(x) abind(x[uids[i]]))
      tmp[[i]]<-do.call("gabind",tmp[[i]])
   }
  names(tmp)<-uids
  return(tmp)
}
