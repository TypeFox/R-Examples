
preproc.crreg <- function(  data = data  , m = m  , trans = trans , status = status , cens.code = cens.code ){

  rownames(data)<-seq_len(nrow(data))
  I <- data[,'left'] != data[,'right' ] & data[,'right'] != Inf & data[,status] != cens.code
  data2<-data[I,]
  dim(data2)
  dataE<-data[!rownames(data)%in%rownames(data2),]
  or<-order(c(as.numeric(rownames(data2)),as.numeric(rownames(dataE))))
  data1<-t(apply( dataE , 1 , function(x) as.numeric(rep(x['left'] , m ))))
dim(data1)[1]+dim(data2)[1]
return( list( data2 = data2 , data1 = data1 , or = or, I = I ) )
}
