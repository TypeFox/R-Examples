gen.sformlist<-function(evl,sforms,cond=FALSE,interval=FALSE,warn=TRUE,...){
   inform<-names(evl$eventlist)
   sformlist<-list()
   if(warn){if(any(grepl("[[:punct:]]|[[:digit:]]",sforms))){if(cond){stop("Sorry, special regex not supported for conditional s-forms. \n \t See ?gen.sformlist for an example of how to do this manually.")};sfin<-sforms[which(grepl("[[:punct:]]|[[:digit:]]",sforms))];cat(paste("Note: \n \t ", sfin, " S-form(s) contains special regex syntax. \n  \t Using slow search methods. \n",sep=""))}}
   for(i in 1:length(inform)){
      evl1<-attr(evl$eventlist[[i]],"char")
      attr(evl1,"a.uid")<-evl$event.key
      if(!is.null(evl$null.events)){
      attr(evl1,"a.null")<-evl$null.events}
      sformlist[[i]]<-mapply(function(x) gen.sform(evl1,x,...),sforms)
   }
  outpm<-sfl2sl(sformlist)
   if(cond){
      suffmat<-lapply(evl$eventlist,function(mx) sapply(sforms,function(my) unlist(regmat.ind(my,attr(mx,"char"))[,2])))
      cond.names<-function(sforms){
      p1.se<-substr(sforms,1,nchar(sforms)-1)
      p2.se<-substr(sforms,nchar(sforms),nchar(sforms))
      paste("(",p1.se,")",p2.se,sep="")
      }
      sforms2<-sapply(sforms,function(x) substr(x,1,nchar(x)-1))
      outpm2<-gen.sformlist(evl,sforms2,warn=FALSE)
      if(interval){
         sforms3<-sapply(sforms2,function(x) substr(x,1,nchar(x)-1))
         outpm3<-gen.sformlist(evl,sforms3,warn=FALSE)
         for(i in 1:length(outpm)){
            for(j in 1:length(sforms)){
               outpm[[i]][,,j]<-abs(outpm[[i]][,,j]-outpm2[[i]][,,j])
               outpm2[[i]][,,j]<-abs(outpm2[[i]][,,j]-outpm3[[i]][,,j])
               #The following conditional updates the last event after absolute differencing
               if(!is.null(unlist(regmat.ind(sforms[j],attr(evl$eventlist[[i]],"char"))[,2]))){
                  ev.up<-unlist(regmat.ind(sforms[j],attr(evl$eventlist[[i]],"char"))[,2])
                  st.up<-attr(evl$eventlist[[i]],"char")[ev.up]
                  outpm[[i]][ev.up,st.up,j]<-1
               }
            }
         outpm[[i]]<-abind(outpm2[[i]],outpm[[i]],new.names=list(NULL,NULL,c(cond.names(sforms2),cond.names(sforms))))
         }        
      }
      if(!interval){
        for(i in 1:length(outpm)){
            for(j in 1:length(sforms)){
               outpm[[i]][,,j]<-abs(outpm[[i]][,,j]-outpm2[[i]][,,j])
               #The following conditional updates the last event after absolute differencing
               if(!is.null(unlist(regmat.ind(sforms[j],attr(evl$eventlist[[i]],"char"))[,2]))){
                  ev.up<-unlist(regmat.ind(sforms[j],attr(evl$eventlist[[i]],"char"))[,2])
                  st.up<-attr(evl$eventlist[[i]],"char")[ev.up]
                  outpm[[i]][ev.up,st.up,j]<-1
               }
            }
         dimnames(outpm[[i]])[[3]]<-cond.names(sforms)         
         } 
      }
   }
   names(outpm)<-names(evl$eventlist)
   return(outpm)
}
