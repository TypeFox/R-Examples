"loadvector" <-
structure(function(filename){
  vv<-read.table(file=filename,skip=1)
  if (dim(vv)[2]==1)
    vv<-vv[[1]]
  vv
}
, comment = "Load vector(s) from file that was produced by Pajek")
