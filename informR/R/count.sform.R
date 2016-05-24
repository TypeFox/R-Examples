count.sform<-function(evls,sform,new.name=TRUE){
   evls2<-lapply(evls$eventlist,attr,which="char")
   b.str.lst<-lapply(evls2,function(x) regmat.ind(sform,x))
   if(new.name){sfn<-sf2nms(evls$event.key,sform)} else{sfn<-sform}
   outpm<-lapply(b.str.lst,nrow)
   outpm[which(is.null(unlist(outpm)))]<-0
   print(paste("Total instances of",sfn,":",sum(unlist(outpm))),quote=FALSE)
   invisible(outpm)
}
