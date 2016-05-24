quantMask<-function(x,quant=.7){
    if(length(dim(x))>3){x<-colMedian(arrMat(x))}
    if(length(dim(x))==2){x<-colMedian(x)}
    if(length(dim(x))==3){x<-as.vector(x)}
    if(is.vector(x)){maskThresh<-quantile(x,quant)}
        maskvec<-ifelse(x <= maskThresh,0,1)
    return(maskvec)
  }

