random.subset <-
function(F_,L_,gamma, persistence=1000, minimum.class.size=2,replace){
#random.subset##random.subset##random.subset##random.subset##random.subset##random.subset##random.subset##random.subset##random.subset##random.subset#
	### A random subset of the samples is selected.
	#____________________________________________________________________________________________________________________________________
	# INPUT: 	F is the features matrix.
	#			L should be the vector of lebels ordered according to rows of F. The lables are -1 or +1 for negative and positive instances.
	#			persistence:	maximum.num.of.tries.for_randomly.choosing.samples (=1000 as default),
	#							If we try this many times and the optained lables are all the same,
	#							we give up (maybe the WHOLE lables are the same) with an error message.
	#			gamma:	=3/4 by default. gamma fraction of samples are considered as training ones.
	#____________________________________________________________________________________________________________________________________

	# Reading input
	total.samples <- rownames(F_)
	remainder.samples <- c()
	### randomly choosing some samples:
	number.of.tries <- 0
	repeat{ #until 						# This subset contains enough samples from both negative and positive ones.
										# Otherwise, fitting a modle dose not make sense.
										# Of course, we will not persist for ever.
		number.of.tries <- number.of.tries +1
												
		#pick some random samples
		selected.samples.indices <- sample(1:length(total.samples), gamma*length(total.samples),replace=replace)
			# We need to look at the indices because the names might be repeated (e.g. by balancing step).
		selected.samples <- total.samples[selected.samples.indices]
		remainder.samples.indices <- which(is.na(match(1:length(total.samples),selected.samples.indices))) #the rest
		remainder.samples <- total.samples[remainder.samples.indices] 
		# updating features matrix
		X_ <- ignore.redundant(F_[selected.samples, ], num.of.values=1)
		Y_ <- L_[rownames(X_)]
			# Y_ needs to be updated to represent the labels for this set of samples 

		positive.class.size <- length(which(Y_ >= 1))
		negative.class.size <- length(which(Y_ < 1))
		if(min(negative.class.size,positive.class.size) >= minimum.class.size )	
			break
				# For instance, for minimum.class.size=1, 
				# it means at least two samples have different labels => enough. 
		if(number.of.tries >= persistence){
			stop(paste("Not enough variation in the labels!",
						"\n We tried: ", number.of.tries,"times.",
						"\n You may want to look at data more carefully or/to increase:",
						" persistence.") )
		}#End if.												
	}#End repeat.

	return(list(X_=X_, Y_=Y_, remainder.samples=remainder.samples))
#random.subset##random.subset##random.subset##random.subset##random.subset##random.subset##random.subset##random.subset##random.subset##random.subset#	
}#End random.subset <- function.

