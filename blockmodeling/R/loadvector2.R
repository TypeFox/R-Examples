"loadvector2" <-
structure(function(filename){
  if(is.character(filename)) {file<-file(description=filename,open="r")
  }else file<-filename
  nn <-read.table(file=file,nrows=1)
  vv<-read.table(file=file,nrows=nn[[2]])
  if (dim(vv)[2]==1)
    vv<-vv[[1]]
  vv
}
, comment = "Load vector(s) from file that was produced by Pajek")
