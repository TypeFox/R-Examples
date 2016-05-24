scale_MinMax <-
function(x){ 
  tmp_1<-dim(x)[2]
  x_scale<-matrix(0,dim(x)[1],dim(x)[2])
  for(i in 1:tmp_1){
    x_scale[,i]<-(x[,i]-min(x[,i]))/(max(x[,i]-min(x[,i])))
  }
  colnames(x_scale)<-colnames(x)
  rownames(x_scale)<-rownames(x)
  return
  x_scale
}
