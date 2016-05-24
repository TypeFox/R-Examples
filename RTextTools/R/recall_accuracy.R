recall_accuracy <- function(true_labels, predicted_labels) 
{
	true_labels <- as.vector(true_labels)
	predicted_labels <- as.vector(predicted_labels,mode=class(true_labels))
	analyze <- predicted_labels == true_labels

	accuracy <- length(analyze[analyze == TRUE])/length(true_labels)
	return(accuracy)
}