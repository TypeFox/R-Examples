order.fdata<-function(y,fdataobj,na.last = TRUE, decreasing = FALSE){
 or<-order(y,na.last=na.last,decreasing=decreasing)
 fdataobj$data<-fdataobj$data[or,]
 fdataobj
}

