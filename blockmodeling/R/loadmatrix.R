"loadmatrix" <-
structure(function(filename){
  if(is.character(filename)) {file<-file(description=filename,open="r")
  }else file<-filename
nn<-read.table(file=file,nrows=1)
if (length(nn) == 2)
  { xx<-read.table(file=file,nrows=nn[[2]],fill=TRUE)
    n<-read.table(file=file,skip=1,nrows=nn[[2]])
    n<-as.matrix(n)
    rownames(n)<-xx[[2]]
    colnames(n)<-xx[[2]] }
 else
   {xxrow<-read.table(file=file,nrows=nn[[3]],fill=TRUE)
    xxcol<-read.table(file=file,nrows=nn[[2]]-nn[[3]],fill=TRUE)
    n<-read.table(file=file,skip=1,nrows=nn[[3]])
    n<-as.matrix(n)
    rownames(n)<-xxrow[[2]]
    colnames(n)<-xxcol[[2]] }
  as.matrix(n)
  }
, comment = "Load matrix from file that was produced by Pajek")
