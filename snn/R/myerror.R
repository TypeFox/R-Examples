myerror <-
function(predict, true){
	# compute the error of the predict list given the true list
	
	if(length(predict) != length(true)) stop("predict and true have different length")	
	ntest = length(predict)	
	error = length(which(predict != true))/ntest
	error
	
}
