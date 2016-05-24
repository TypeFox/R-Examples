write.sym.table <-
function(sym.data,file,sep,dec,row.names=NULL,col.names=NULL) {
  write.table(sym.data$meta,file,sep=as.character(sep),dec=dec,quote=FALSE,
              row.names=c(row.names),col.names=c(col.names))
}
