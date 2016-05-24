alpha2h <-
function(alpha,  at, xseq){
  if(alpha>1|alpha<0) stop("Specified alpha should be  between 0 and 1.")
  n<-length(xseq)
  nn<-round(alpha*n)
  ds<-abs(xseq-at)  
  h<-sort(ds, decreasing=T)[nn]
  return(h) 
 }
