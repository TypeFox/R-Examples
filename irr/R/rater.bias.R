# rater.bias computes a Chi-squared value for a systematic bias
# of one rater compared with another.

rater.bias<-function(x) {
 if(is.matrix(x) || is.data.frame(x)) {
  dimx<-dim(x)
  if(length(dimx) == 2) {
   # if dimension lengths are unequal, assume it's a nx2 score matrix
   if(dimx[1] != dimx[2]) {
    if(dimx[1] == 2) rbx<-as.matrix(table(x[1,],x[2,]))
    else rbx<-as.matrix(table(x[,1],x[,2]))
   }
   else rbx<-as.matrix(x)
   # print(rbx)
   rbb<-sum(rbx[upper.tri(rbx)])
   rbc<-sum(rbx[lower.tri(rbx)])
   rb <-abs(rbb/(rbb+rbc))
   rbstat<-(rbb-rbc)^2/(rbb+rbc)
   rblist<-structure(list(method="Rater bias coefficient",subjects=sum(rbx),
    raters=2,irr.name="Ratio",value=rb,stat.name="Chisq(1)",
    statistic=rbstat,p.value=1-pchisq(rbstat,1)),class="irrlist")
   return(rblist)
  }
  else cat("Dimension higher than 2, cannot compute\n")  
 }
 cat("Usage: rater.bias(x)\n")
 cat("\twhere x is an nx2 or 2xn  matrix of category scores for n objects\n")
 cat("\tor a CxC matrix or data frame of rater agreement on C categories\n")
}
