generate.sym.table <-
function(sym.data,file,sep,dec,row.names=NULL,col.names=NULL) {
  ncolumn<-rep(sym.data$sym.var.types[1],sym.data$N)
  temp.data <- cbind(temp=ncolumn,sym.data$data[,1:ncol(sym.data$data)])
  colnames(temp.data)[1]<-sym.data$sym.var.types[1]
  pos<-1 # pos = valid number of columns inserted
  if((sym.data$sym.var.types[1]=='$C')||(sym.data$sym.var.types[1]=='$c')) {
    pos<-2
  }
  if((sym.data$sym.var.types[1]=='$I')||(sym.data$sym.var.types[1]=='$i')) {
    pos<-3
  }
  if((sym.data$sym.var.types[1]=='$H')||(sym.data$sym.var.types[1]=='$h')||
       (sym.data$sym.var.types[1]=='$S')||(sym.data$sym.var.types[1]=='$s')) {
    ncolumn<-rep(sym.data$sym.var.length[1],sym.data$N)
    temp.data<-cbind(temp.data[,1:pos],temp=ncolumn,temp.data[,(pos+1):ncol(temp.data)])
    colnames(temp.data)[(pos+1)]<-sym.data$sym.var.names[1]
    pos<-pos+1+sym.data$sym.var.length[1]
  }    
  for(j in 2:sym.data$M)  {
    ncolumn<-rep(sym.data$sym.var.types[j],sym.data$N)
    temp.data<-cbind(temp.data[,1:pos],temp=ncolumn,temp.data[,(pos+1):ncol(temp.data)])
    colnames(temp.data)[pos+1]<-sym.data$sym.var.types[j]
    pos<-pos+1
    if((sym.data$sym.var.types[j]=='$C')||(sym.data$sym.var.types[j]=='$c')) {
      pos<-pos+1
    }
    if((sym.data$sym.var.types[j]=='$I')||(sym.data$sym.var.types[j]=='$i')) {
      pos<-pos+2
    }    
    if((sym.data$sym.var.types[j]=='$H')||(sym.data$sym.var.types[j]=='$h')||
         (sym.data$sym.var.types[j]=='$S')||(sym.data$sym.var.types[j]=='$s')) {
      ncolumn<-rep(sym.data$sym.var.length[j],sym.data$N)
      temp.data<-cbind(temp.data[,1:pos],temp=ncolumn,temp.data[,(pos+1):ncol(temp.data)])
      colnames(temp.data)[(pos+1)]<-sym.data$sym.var.names[j]
      pos<-pos+1+sym.data$sym.var.length[j]
    }
  }
  write.table(temp.data,file,sep=as.character(sep),dec=dec,quote=FALSE,
              row.names=c(row.names),col.names=c(col.names))
}
