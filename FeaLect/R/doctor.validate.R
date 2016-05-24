doctor.validate <-
function(true.labels,predictions){
### validating: vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# validates the trained doctor as a model and compares the errors.
	
	# Input check:
	if(length(true.labels)!=length(predictions))	
		stop("predictions should be as long as true labels.")
	if(prod(as.numeric( names(true.labels)==names(predictions))) ==0 )	
		stop("predictions and true labels should have been named in the same order.")
		
	# These are indices, not names.
	positives <- which(true.labels > 0)
	negatives <- which(true.labels <= 0)
	predicted.positives <- which(predictions > 0)
	predicted.negatives <- which(predictions <= 0)
	true.positives <- intersect(positives,predicted.positives)
	true.negatives <- intersect(negatives,predicted.negatives)
	precision <- length(true.positives)/(length(predicted.positives))
	if(is.na(precision))	
		precision <- 0	# because some devide by zero has happend.				
	recall <- length(true.positives)/(length(positives))				
	if(is.na(recall))	
		recall <- 0	# because some devide by zero has happend.
	f.measure <- 2*precision*recall/(precision+recall)
	if(is.na(f.measure))	
		f.measure <- 0	# because some devide by zero has happend.
			
	#Keeping track of the miss-labeled samples,
	false.positives <- names(predictions)[intersect(predicted.positives, negatives)]
	false.negatives <- names(predictions)[intersect(predicted.negatives, positives)]
	mislabeled <- union(false.positives,false.negatives)

	# Result:
	result <- list(f.measure=f.measure,mislabeled=mislabeled,precision=precision,recall=recall)
		
	return(result)
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv	
}#End validate <- function.	

