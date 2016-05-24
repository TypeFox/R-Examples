#calculate a randomization-based interval for the test statistic based on empirical randomization distribution
#input: results object outputted by randTest(), coverage level (defaults to .95),
#output: randomization-based interval for the test statistic or vector of test statistics

randInterval=function(results,coverage=.95){
	if(length(results$obs_stat)==1){
		if(results$alternative=="two.sided"){
			c(quantile(results$perm_stats,(1-coverage)/2),quantile(results$perm_stats,(1+coverage)/2))
		}else if(results$alternative=="greater"){
			c(quantile(results$perm_stats,0),quantile(results$perm_stats,coverage))
		}else if(results$alternative=="less"){
			c(quantile(results$perm_stats,1-coverage),quantile(results$perm_stats,1))
		}
	}else{
		if(results$alternative=="two.sided"){
			sapply(1:length(results$obs_stat),function(j) c(quantile(results$perm_stats[,j],(1-coverage)/2),quantile(results$perm_stats[,j],(1+coverage)/2)))
		}else if(results$alternative=="greater"){
			sapply(1:length(results$obs_stat),function(j) c(quantile(results$perm_stats[,j],0),quantile(results$perm_stats[,j],coverage)))
		}else if(results$alternative=="less"){
			sapply(1:length(results$obs_stat),function(j) c(quantile(results$perm_stats[,j],1-coverage),quantile(results$perm_stats[,j],1)))
		}
	}
}
