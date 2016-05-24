measure.edit<-function(file,plate){
  if(file=="") file<-paste(".last",plate,"_measure.txt",sep="")
  measure<-read.delim(file,header=FALSE)
  for(i in 1:length(measure[1,])) measure[,i]<-as.character(measure[,i])
  names(measure)[1]<-"Compound"
  for(con in 1:(length(measure[1,])-1)) names(measure)[con+1]<-paste("Con",con,sep="")
  measure<-edit(measure)
  write.table(measure,file=file,sep="\t",row.names=FALSE,col.names=FALSE)
}

control.edit<-function(file,plate){
  if(file=="") file<-paste(".last",plate,"_control.txt",sep="")
  control<-read.delim(file,header=FALSE)
  for(i in 1:length(control[1,])) control[,i]<-as.character(control[,i])
  names(control)[1]<-"Compound"
  for(con in 1:(length(control[1,])-1)) names(control)[con+1]<-paste("Con",con,sep="")
  control<-edit(control)
  write.table(control,file=file,sep="\t",row.names=FALSE,col.names=FALSE)  
}

dilution.edit<-function(file,plate){
  if(file=="") file<-paste(".last",plate,"_dilution.txt",sep="")
  dilution<-read.delim(file,header=FALSE)
  for(i in 1:length(dilution[1,])) dilution[,i]<-as.character(dilution[,i])
  names(dilution)[1]<-"Compound"
  for(con in 1:(length(dilution[1,])-1)) names(dilution)[con+1]<-paste("Con",con,sep="")
  dilution<-edit(dilution)
  write.table(dilution,file=file,sep="\t",row.names=FALSE,col.names=FALSE)
}
