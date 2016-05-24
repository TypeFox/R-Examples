`slbind`<-
function (sformstats, statslist, type = 1, new.names=FALSE,...) 
{
    if (length(statslist) != length(sformstats)) {
        stop("\n Both objects must be the same length.")
    }
    newstatsl <- statslist
    for (i in 1:length(statslist)) {
        newstatsl[[i]][[type]] <- abind(statslist[[i]][[type]], 
            sformstats[[i]])
    }
    if(is.logical(new.names)){
       if(!new.names){return(newstatsl)}
       if(new.names){
          for(i in 1:length(newstatsl)){
          b1<-dimnames(statslist[[i]][[type]])[[3]]
          b2<-dimnames(newstatsl[[i]][[type]])[[3]]
          b3<-b2[which(!b2%in%b1)]
          dimnames(newstatsl[[i]][[type]])[[3]]<-c(b1,sapply(b3,sf2nms,...))
         }
        return(newstatsl)
      }
    }
    if(is.character(new.names)){
       for(i in 1:length(newstatsl)){
       b1<-dimnames(statslist[[i]][[type]])[[3]]
       dimnames(newstatsl[[i]][[type]])[[3]]<-c(b1,new.names)
      }
     return(newstatsl)
    }
}
