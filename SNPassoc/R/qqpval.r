qqpval<-function(p, pch=16, col=4, ...){
   p<-p[!is.na(p)]
   n<-length(p)
   pexp<-(1:n)/(n+1)
   plot(-log(pexp,10), -log(sort(p),10), xlab="-log(expected P value)", ylab="-log(observed P value)", pch=pch, col=col, ...)
   abline(0,1,col=2)
}
