symscorestat<-function(y,scores=NULL,exact=F,sides=1){
   if(is.null(scores)) scores<-seq(length(y))
   scores<-scores[order(abs(y))]
   if(exact){
      out<-.Fortran("signtestperm",as.double(y),as.double(scores),as.integer(length(y)),out=as.integer(0),PACKAGE="MultNonParam")$out
#     cat("out",out)
      pvo<-out*2^(-length(y))
   }else{
      zstat<- -(sum(scores[y>0])-sum(scores)/2)/
         sqrt(sum(scores^2)/4)
      pvo<-pnorm(-zstat)
   }
   if(sides==2) pvo<-2*min(pvo,1-pvo)
   return(list(pv=pvo,tot=2^length(y)))
}
