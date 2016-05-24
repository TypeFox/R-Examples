#plot observed test statistic against empirical randomization distribution
#input: results object outputted by randTest(), coverage level (defaults to .95),
#	breaks in the histogram (breaks), plot layout dimensions if multiple plots (plotDim)
#output: histogram(s) of permuted test statistics,
#	observed test statistic(s) demarcated by red vertical line,
#	randomization-based interval(s) demarcated by blue vertical lines

randPlot=function(results,coverage=.95,breaks=10,plotDim=c(length(results$obs_stat),1)){
	interval=randInterval(results,coverage)
	if(length(results$obs_stat)==1){
		par(mfrow=plotDim)
		hist(results$perm_stats,breaks=breaks,xlab="Test Statistic",main="Randomization Distribution of Test Statistic")
		abline(v=results$obs_stat,col="red",lwd=5)
		abline(v=results$null,col="black",lwd=3)
		abline(v=c(interval[1],interval[2]),col="blue",lwd=3,lty=2)
	}else{
		par(mfrow=plotDim)		
		for(j in 1:length(results$obs_stat)){
			hist(results$perm_stats[,j],breaks=breaks,xlab="Test Statistic",main=paste("Randomization Distribution of Test Statistic",j))
			abline(v=results$obs_stat[j],col="red",lwd=5)
			abline(v=results$null[j],col="black",lwd=3)
			abline(v=c(interval[1,j],interval[2,j]),col="blue",lwd=3,lty=2)
		}
	}
	par(mfrow=c(1,1))
}
