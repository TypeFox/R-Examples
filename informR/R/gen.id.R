`gen.id` <-
function(rawevents,print=TRUE){
   ea<-rawevents
   ea.u<-unique(ea)
   if(length(ea.u)>52){stop("More than 52 event classes currently not supported")}
   stemids<-c(letters,LETTERS)
   ea.id<-cbind(id=stemids[1:length(ea.u)],event.type=ea.u)
   if(print){print(ea.id)}
   eventlist<-ea.id[match(ea,ea.id[,2]),1]
   attr(eventlist,"event.key")<-ea.id
   eventlist
}
