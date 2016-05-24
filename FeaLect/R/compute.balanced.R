compute.balanced <-
function(F_,L_){
	# balancing:
	# 	If the size of classes are not equal according to L, 
	# 	we balance the (F,L) system by adding (oversampling) or removing (unersampling) some samples.
	#________________________________________________________________________________________________

	# Computing new indices:
	positive.indecies <- which(L_ >=1)
	negative.indecies <- which(L_ < 1)
	size.difference <- length(negative.indecies) - length (positive.indecies)
	if (size.difference  >= 0){	# more negatives,
		small.class <- positive.indecies
		big.class <- negative.indecies
	} else {
		small.class <- negative.indecies
		big.class <- positive.indecies
	}#End elsif.
		
	new.indices <- sample(small.class, abs(size.difference),replace = TRUE)
	# Adding new indices
	F_update <- rbind(F_, F_[new.indices,])
	L_update <- c(L_, L_[new.indices])	
	#message(dim(F_))

	return(list(F_=F_update,L_=L_update))
}#End compute.balanced <- function(F,Y).

