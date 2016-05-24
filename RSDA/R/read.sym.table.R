read.sym.table <-
function(file,header=TRUE,sep,dec,row.names=NULL) {
  meta.data<-read.table(file,header,sep=as.character(sep),dec=as.character(dec),row.names=c(row.names),check.names=FALSE)
  n.sym.objects<-dim(meta.data)[1]
  meta.M<-dim(meta.data)[2]
  sym.var.types<-list()
  sym.var.length<-rep(0,length(meta.M))
  sym.var.names<-list()
  sym.var.starts<-list()
  sym.obj.names<-rownames(meta.data)
  del.columns<-c()
  del.columns.length<-0
  if(header==TRUE) 
    meta.types<-colnames(meta.data)
  else 
    stop("Data file have to have a header")
  meta.types.orig<-meta.types
  for(i in 1:length(meta.types)) {
    meta.types[i]<-substr(meta.types[i],start=1,stop=2) 
  }
  for(j in 1:length(meta.types)) {    
    if((meta.types[j]=='$C')||(meta.types[j]=='$c')) {
      sym.var.types[j]<-'$C'
      sym.var.length[j]<-1
      sym.var.names[j]<-meta.types.orig[j+1]
      sym.var.starts[j]<-j+1
    }   
    else 
      if((meta.types[j]=='$I')||(meta.types[j]=='$i')) {
        sym.var.types[j]<-'$I'
        sym.var.length[j]<-2
        sym.var.names[j]<-meta.types.orig[j+1]
        sym.var.starts[j]<-j+1
      }  
    else 
      if((meta.types[j]=='$H')||(meta.types[j]=='$h')) {
        sym.var.types[j]<-'$H' 
        sym.var.length[j]<-as.integer(meta.data[2,j+1])
        del.columns[del.columns.length+1]<-j+1
        del.columns.length<-del.columns.length+1
        sym.var.names[j]<-meta.types.orig[j+1]
        sym.var.starts[j]<-j+2
      }  
    else
      if((meta.types[j]=='$S')||(meta.types[j]=='$s')) {
        sym.var.types[j]<-'$S'
        sym.var.length[j]<-as.integer(meta.data[2,j+1])
        del.columns[del.columns.length+1]<-j+1
        del.columns.length<-del.columns.length+1 
        sym.var.names[j]<-meta.types.orig[j+1]
        sym.var.starts[j]<-j+2
      }   
      else
        sym.var.types[j]<-'NA'      
  }  
  del1<-match(sym.var.types,c('$C','$I','$H','$S'),0)
  for(k in 1:del.columns.length) {
    sym.var.types[del.columns[k]]<-'$H'
  }
  del2<-match(sym.var.types,c('$C','$I','$H','$S'),0)
  sym.var.types<-sym.var.types[del1>0]
  sym.var.length<-sym.var.length[del1>0]
  n.sym.var<-length(sym.var.length)
  data.matrix<-as.data.frame(meta.data[,del2==0])
  sym.data<-list(N=n.sym.objects,M=n.sym.var,sym.obj.names=sym.obj.names,
                 sym.var.names=unlist(sym.var.names),sym.var.types=unlist(sym.var.types),
                 sym.var.length=sym.var.length,sym.var.starts=unlist(sym.var.starts),
                 meta=meta.data,data=data.matrix)
  return(sym.data)
}
