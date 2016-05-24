`gen.intercepts` <-
function(evl, basecat=NULL, type=1, contr=TRUE){
   contr.evs<-evl$event.key[,2]
   if(!is.null(evl$null.events)){
   contr.evs<-evl$event.key[-which(evl$event.key[,2]%in%evl$null.events),2]
   }
   val.events<-lapply(evl$eventlist,function(x) evl$event.key[,2][match(attr(x,"char"),evl$event.key[,1])])
   if(contr){
   if(is.null(basecat)){basecat<-sort(contr.evs)[1]}
   baselines<-contr.SAS(contr.evs,contrasts=FALSE)
   baselines<-baselines[,-which(colnames(baselines)==basecat)]}
   if(!contr){
   baselines<-contr.SAS(contr.evs,contrasts=FALSE)
   }

   inform<-names(evl$eventlist)
   statsl<-list()
   for(i in 1:length(inform)){
      evl2<-val.events[[i]]
      if(type==1){
         statsl[[i]]<-list(global=array(0,dim=c(length(evl2),dim(baselines)[[1]],dim(baselines)[[2]])))
         for(j in 1:length(evl2)){
            statsl[[i]]$global[j,,]<-baselines
         }
         dimnames(statsl[[i]]$global)<-list(1:length(evl2),dimnames(baselines)[[1]],dimnames(baselines)[[2]])
      }
      if(type==2){
         statsl[[i]]<-list(local=array(0,dim=c(length(evl2),dim(baselines)[[1]],dim(baselines)[[2]])))
         for(j in 1:length(evl2)){
            statsl[[i]]$local[j,,]<-baselines
         }
         dimnames(statsl[[i]]$local)<-list(1:length(evl2),dimnames(baselines)[[1]],dimnames(baselines)[[2]])
      }
   }
 statsl
}
