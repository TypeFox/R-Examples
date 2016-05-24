dA<-function(w,A,dw){
wa<-sqrt(sum((w*(A%*%w))))
dummy<-(1/wa)*(diag(length(w))- w%*%t(w)%*%A/(wa^2))%*%dw
return(dummy)
}