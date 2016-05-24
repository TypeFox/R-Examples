maxent <-
function(feature_matrix,code_vector,l1_regularizer=0.0,l2_regularizer=0.0,use_sgd=FALSE,set_heldout=0,verbose=FALSE) {	
    suppressWarnings(sink())
    feature_matrix <- as.compressed.matrix(feature_matrix);
	train_maxent <- function(feature_matrix,code_vector,l1_regularizer=0.0,l2_regularizer=0.0,use_sgd=FALSE,sgd_iter=30,sgd_eta0=1.0,sgd_alpha=0.85,set_heldout=0) {
        
        if (length(unique(feature_matrix@ja)) > 1) { ja <- sapply(feature_matrix@ja,toString); }
        else { ja <- sapply(feature_matrix@ra,toString); }
		code_vector <- sapply(code_vector,toString);
		maximumentropy$add_samples(as.integer(feature_matrix@dimension[1]),as.integer(feature_matrix@dimension[2]),code_vector,feature_matrix@ia,ja,feature_matrix@ra);
		model <- maximumentropy$train_model(l1_regularizer,l2_regularizer,use_sgd,sgd_iter,sgd_eta0,sgd_alpha,set_heldout);
		
		return(model);
	}
	
	if (l1_regularizer > 0 && l2_regularizer > 0) {
		print("ERROR: L1 and L2 regularization cannot be used together.");
		return(NULL);
	}
	
    if (l2_regularizer > 0 && use_sgd==TRUE) {
        print("ERROR: L2 regularization is currently not supported in SGD mode.");
		return(NULL);
    }
    
    if (length(unique(code_vector)) > 255) {
        print("ERROR: Too many types of labels (>255 unique labels).");
        return(NULL);
    }

    if (verbose == FALSE) {
        if(.Platform$OS.type == "unix") {
            sink("/dev/null")
        } else {
            sink("NUL")
        }
    }
    
    model <- train_maxent(feature_matrix,code_vector,l1_regularizer,l2_regularizer,use_sgd,set_heldout);
	
	weights <- as.data.frame(cbind(model[[2]],model[[3]],model[[4]]));
	colnames(weights) <- c("Weight","Label","Feature");
	container <- new("maxent", model=model[[1]], weights=weights);
    
    suppressWarnings(sink())

	return(container);
}