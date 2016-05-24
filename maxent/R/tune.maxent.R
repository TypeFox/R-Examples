tune.maxent <-
function(feature_matrix,code_vector,nfold=3,showall=FALSE,verbose=FALSE) {
    recall_accuracy <- function(true_labels, predicted_labels) 
    {
        true_labels <- as.vector(true_labels)
        predicted_labels <- as.vector(predicted_labels,mode=class(true_labels))
        analyze <- predicted_labels == true_labels

        accuracy <- length(analyze[analyze == TRUE])/length(true_labels)
        return(accuracy)
    }
    
    feature_matrix <- as.compressed.matrix(feature_matrix)
    if (verbose==TRUE) cat("Testing 18 parameter configurations...\n")
    l1_params <- c(seq(0,1,0.2),rep(0,length(seq(0,1,0.2))),seq(0,1,0.2))
    l2_params <- c(rep(0,length(seq(0,1,0.2))),seq(0,1,0.2),rep(0,length(seq(0,1,0.2))))
    sgd_params <- c(rep(FALSE,12),rep(TRUE,6))
    set_heldout_params <- c(rep(0,6),round(dim(feature_matrix)[1]/(nfold*2)),rep(0,5),rep(0,6))
    
    rand <- sample(nfold,dim(feature_matrix)[1], replace=TRUE)
    fit_accuracy <- c()
    
    for (n in 1:length(l1_params)) {
        cv_accuracy <- c()
        for (i in sort(unique(rand))) {
            model <- maxent(feature_matrix[rand!=i,],code_vector[rand!=i],l1_regularizer=l1_params[n],l2_regularizer=l2_params[n],use_sgd=sgd_params[n],set_heldout=set_heldout_params[n],verbose=FALSE)
            pred <- predict(model,feature_matrix[rand==i,])
            pred <- pred[,1]

            cv_accuracy[i] <- recall_accuracy(code_vector[rand==i],pred)
        }
        
        if (verbose==TRUE) cat(paste("Configuration: ",n,"\tAccuracy (",nfold,"-fold cross-validation): ",mean(cv_accuracy),"\n",sep=""))
        fit_accuracy[n] <- mean(cv_accuracy)
    }
    
    names <- c("l1_regularizer","l2_regularizer","use_sgd","set_heldout","accuracy","pct_best_fit")
    if (showall==FALSE) {
        optimal_fit <- which.max(fit_accuracy)
        values <- c(l1_params[optimal_fit],l2_params[optimal_fit],sgd_params[optimal_fit],set_heldout_params[optimal_fit],max(fit_accuracy),100)
        optimal <- rbind(values,deparse.level=0)
        colnames(optimal) <- names
    } else {
        optimal <- rbind()
        best <- max(fit_accuracy)
        for (optimal_fit in 1:length(l1_params)) {
            values <- c(l1_params[optimal_fit],l2_params[optimal_fit],sgd_params[optimal_fit],set_heldout_params[optimal_fit],fit_accuracy[optimal_fit],fit_accuracy[optimal_fit]/best)
            optimal <- rbind(optimal,values,deparse.level=0)
        }
        colnames(optimal) <- names
    }
    
    return(optimal)
}