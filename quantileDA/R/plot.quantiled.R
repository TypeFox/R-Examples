plot.quantiled<-function(x,...)
{
	out=x
	nomi<-c("quantile classifier","median classifier","centroid classifier")
	nfold<-length(out$folds)
	y<-as.factor(rep(nomi,each=nfold))
	x<-c(out$me.test,out$me.median,out$me.centroid)
	plot(x~y,col=2,xlab="",ylab="CV misclassification errors")
}