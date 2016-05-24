decimal.align<-function(x,dechar=".",nint=NA,ndec=NA,pad.left=TRUE) {

 splitchar<-paste("[",dechar,"]",sep="")
 splitlist<-strsplit(as.character(x),splitchar)
 ints<-unlist(lapply(splitlist,"[",1))
 ints[is.na(ints)]<-"0"
 if(pad.left) {
  if(is.na(nint)) nint<-max(nzchar(ints))
  ints<-sprintf(paste("%",nint,"s",sep=""),ints)
 }
 if(is.na(ndec)) ndec<-max(nzchar(unlist(lapply(splitlist,"[",2))))
 decs<-unlist(lapply(splitlist,"[",2))
 decs[is.na(decs)]<-"0"
 decs<-sprintf(paste("%-",ndec,"s",sep=""),
  decs)
 return(paste(ints,decs,sep=dechar))
}
