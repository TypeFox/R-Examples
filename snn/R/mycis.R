mycis <-
function(predict1, predict2){
	# compute the classification instability (cis) of a classification procedure given its two predicted lists trained from two i.i.d. training data sets.
	
	if(length(predict1) != length(predict2)) stop("two predict lists have different length")	
	ntest = length(predict1)	
	cis = length(which(predict1 != predict2))/ntest
	cis

}
