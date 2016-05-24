`hsWriteTable` <-
function(d, file="",sep='\t') {
  write.table(file=file,d,quote=FALSE,sep=sep,row.names=FALSE,col.names=FALSE)
}
