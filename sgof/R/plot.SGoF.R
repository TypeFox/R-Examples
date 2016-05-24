plot.SGoF <-
function(x, ...){
data <- x$data
Adjusted.pvalues<-x$Adjusted.pvalues
plot(sort(data),sort(Adjusted.pvalues),main="SGoF Adjusted p-values",xlim=c(0,1),ylim=c(0,1),xlab="Unadjusted p-values",type="p",ylab="Adjusted p-values",...)
abline(0,1,lty=2,lwd=1.5)
}
