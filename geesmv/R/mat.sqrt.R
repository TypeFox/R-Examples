mat.sqrt <-
function(A)
{
   ei<-eigen(A)
   d<-ei$values
   d<-(d+abs(d))/2
   d2<-sqrt(d)
   ans<-ei$vectors %*% diag(d2) %*% t(ei$vectors)
   return(ans)
}
