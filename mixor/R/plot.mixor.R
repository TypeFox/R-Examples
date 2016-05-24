plot.mixor <-
function(x, ...) {
	par(mfrow=c(dim(x$EBmean)[2],2))
	for (i in 1:dim(x$EBmean)[2]) {
		hist(x$EBmean[,i], xlab="Empirical Bayes Estimates", ylab="Frequency", main=paste("EBmean",i,sep=" "))
		qqnorm(x$EBmean[,i])
		qqline(x$EBmean[,i])
   	}
}