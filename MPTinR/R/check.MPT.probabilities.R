
.check.MPT.probabilities <- function(tree){
	tmp.env <- new.env()
	temp.param.names <- .find.MPT.params(tree)
	temp.param.val <- runif(length(temp.param.names))
	temp.branch <- sapply(tree,length)
	prob <- rep(NA,length(temp.branch))
   
	for (i in 1:length(temp.param.val)) {
		assign(temp.param.names[i],temp.param.val[i], envir = tmp.env)
	}
	temp.check <- sapply(unlist(tree),eval, envir = tmp.env)
	for (i in 1:length(temp.branch)){
		if (i==1) prob[1] <- sum(temp.check[1:temp.branch[1]])
		else prob[i] <- sum(temp.check[(sum(temp.branch[1:(i-1)])+1):sum(temp.branch[1:i])])
	}
	prob <- round(prob,digits=6)
	return(prob)
}
