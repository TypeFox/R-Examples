mat.sqrt.inv <-
function(A)
{
   ei<-eigen(A)
   d<-ei$values
   d<-(d+abs(d))/2
   d2<-1 / sqrt(d)
   d2[d == 0]<-0
   ans<-ei$vectors %*% diag(d2) %*% t(ei$vectors)
   return(ans)
}
