recluster.hist <-function(x) {
	perc_repl<-100*(length(x)-length(unique(x)))/length(x)
	zero<-sum(x==0)
	hist(x,breaks=c(-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),main=paste("Zero x=",zero," Percentage of tied cells=", round(perc_repl,2)))
}
