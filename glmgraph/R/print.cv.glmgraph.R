print.cv.glmgraph <- function(x,...){
	cat("The solution path for the regularization parameter lambda1 and lambda2 is:\n")
	print(x$cvmat)
	if(x$type.measure=="mse" || x$type.measure=="mae" || x$type.measure=="deviance"){
		cat(paste("The minimum ",x$type.measure,"value",x$cvmin," is achieved when lambda2 is: ",x$lambda2.min," and lambda1 is: ",x$lambda1.min,"\n"))
	}else if(x$type.measure=="auc"){
		cat(paste("The maximum ",x$type.measure,"value",-x$cvmin," is achieved when lambda2 is: ",x$lambda2.min," and lambda1 is: ",x$lambda1.min,"\n"))
	}
}
