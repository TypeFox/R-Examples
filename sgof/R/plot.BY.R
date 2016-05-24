plot.BY <-
function(x, ...){
data <- x$data
Adjusted.pvalues<-x$Adjusted.pvalues
plot(sort(data),sort(Adjusted.pvalues),main="BY Adjusted p-values",ylim=c(0,1),xlim=c(0,1),xlab="Unadjusted p-values",type="p",ylab="Adjusted p-values",...)
abline(0,1,lty=2,lwd=1.5)
}
