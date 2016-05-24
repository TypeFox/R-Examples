`gen.evl` <-
function(eventlist,null.events=NULL){
   el.dim<-NCOL(eventlist)
   if(el.dim==1||el.dim>3){stop(paste("","incorrect number of dimensions in eventlist",sep="\n\t"))}
   if(el.dim<3){
      val.events<-eventlist
      sym.events<-gen.id(val.events[,1],print=FALSE)
      evls<-list()
      evns<-unique(val.events[,1])
      evls.char<-glapply(sym.events,val.events[,2],function(x) x,regroup=FALSE)
      evls<-glapply(val.events[,1],val.events[,2],function(x) match(x,evns),regroup=FALSE)
      if(!is.null(null.events)){
      bout<-match(null.events,evns)
      b1<-attr(sym.events,"event.key")[which(!attr(sym.events,"event.key")[,2]%in%null.events),1]
      evls<-lapply(evls.char,function(x) match(x,b1))
      evls<-lapply(evls,function(x) ifelse(is.na(x),0,x))
      #evls<-lapply(evls,function(x){x[which(x%in%bout)]<-0; x})
      }
     for(i in 1:length(evls)){attr(evls[[i]],"char")<-evls.char[[i]]}
     out.evl<-list()
     out.evl$eventlist<-evls
     out.evl$event.key<-attr(sym.events,"event.key")
     if(!is.null(null.events)){out.evl$null.events<-null.events}
     #Check for event lists will less than two events and warn.
     if(any(lapply(out.evl$eventlist,length)<2)){
        warning(paste("Event History(ies):",which(lapply(out.evl$eventlist,length)<2),"has(have) less than 2 events. Removing from eventlist.",sep=" "))
     out.evl$eventlist[which(lapply(out.evl$eventlist,length)<2)]<-NULL
     }
     return(out.evl)
   }
   if(el.dim==3){
      val.events<-eventlist
      sym.events<-gen.id(val.events[,1],print=FALSE)
      evls<-list()
      evns<-unique(val.events[,1])
      evls.char<-glapply(sym.events,val.events[,3],function(x) x,regroup=FALSE)
      evls<-glapply(val.events[,1],val.events[,3],function(x) match(x,evns),regroup=FALSE)
      evls.time<-glapply(val.events[,2],val.events[,3],function(x) as.numeric(x),regroup=FALSE)
      if(!is.null(null.events)){
      bout<-match(null.events,evns)
      b1<-attr(sym.events,"event.key")[which(!attr(sym.events,"event.key")[,2]%in%null.events),1]
      evls<-lapply(evls.char,function(x) match(x,b1))
      evls<-lapply(evls,function(x) ifelse(is.na(x),0,x))
      #evls<-lapply(evls,function(x){x[which(x%in%bout)]<-0; x})
      }
     for(i in 1:length(evls)){
        evls[[i]]<-cbind(evls[[i]],evls.time[[i]])
        attr(evls[[i]],"char")<-evls.char[[i]]
     }
     out.evl<-list()
     out.evl$eventlist<-evls
     out.evl$event.key<-attr(sym.events,"event.key")
     if(!is.null(null.events)){out.evl$null.events<-null.events}
     #Check for event lists will less than two events and warn.
     if(any(lapply(out.evl$eventlist,nrow)<2)){
        warning(paste("Event History(ies):",which(lapply(out.evl$eventlist,nrow)<2),"has(have) less than 2 events. Removing from eventlist.",sep=" "))
     out.evl$eventlist[which(lapply(out.evl$eventlist,nrow)<2)]<-NULL
     }
     return(out.evl)     
   }
}
