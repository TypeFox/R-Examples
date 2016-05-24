Append <-
function(){
  ssdA<-read.table(file.choose(),header=TRUE,sep=',')
  ssdB<-read.table(file.choose(),header=TRUE,sep=',')
  ssdC<-rbind(ssdA,ssdB)

write.csv(ssdC,file = tclvalue(tcl("tk_getSaveFile")),row.names=FALSE)

}
