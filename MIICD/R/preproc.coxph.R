
preproc.coxph <-
function(  data = data  , k = k ){

  rownames(data)<-seq_len(nrow(data))
  I <-  data[,'left'] != data[,'right'] & data$right != Inf
  data2<-data[I,]
  dataE<-data[!rownames(data)%in%rownames(data2),]
  or<-order(c(as.numeric(rownames(data2)),as.numeric(rownames(dataE))))
  data1<-t(apply( dataE , 1 , function(x) as.numeric(rep(x['left'],k))))
  return(list(data2=data2,data1=data1,or=or,I=I))
}