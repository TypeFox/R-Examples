print.quantileDA<-function(x,...)
{
	out=x
	message("")
	message("Results of the quantile classifier")
	message("")
	cat(paste("Minimum misclassification rates attained at theta = ",out$theta.choice,sep=""),"\n")
	cat(paste("Miclassification rate in the training set: ", round(out$me.train,2),sep=""),"\n")
	cat(paste("Miclassification rate in the test set: ",round(out$me.test,2),sep=""),"\n")
	
}