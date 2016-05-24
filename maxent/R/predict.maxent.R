predict.maxent <-
function(object, feature_matrix, ...) {
    suppressWarnings(sink())
    feature_matrix <- as.compressed.matrix(feature_matrix);
	classify_maxent <-
	function(feature_matrix,model) {
        if (length(unique(feature_matrix@ja)) > 1) { ja <- sapply(feature_matrix@ja,toString); }
        else { ja <- sapply(feature_matrix@ra,toString); }
        results <- maximumentropy$classify_samples(as.integer(feature_matrix@dimension[1]),as.integer(feature_matrix@dimension[2]),feature_matrix@ia,ja,feature_matrix@ra,model);
		
		labels <- as.vector(results[[1]]);
		probabilities <- as.matrix(results[[2]],mode="numeric");
		colnames(probabilities) <- as.vector(results[[3]],mode="character");
		
		return(cbind(labels,probabilities));
	}	

    results <- classify_maxent(feature_matrix,object@model);
	
	return(results);
}