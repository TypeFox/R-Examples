classify_model <-
function(container, model, s=0.01, ...) {

	gc()
    extract_maximum_prob <- function(x) return(x[which.max(x)])
    extract_label_from_prob <- function(x) return(which.max(x))
	extract_label_from_prob_names <- function(x) return(rownames(as.matrix(which.max(x))))
    
    if (pmatch("svm",class(model),nomatch=0) > 0){
        svm_results <- predict(model,container@classification_matrix, prob=TRUE, ...) #Extract Label
        svm_pred <- svm_results[1:length(svm_results)]
        svm_prob <- apply(attr(svm_results,"prob"),1,extract_maximum_prob)

        results_table <- data.frame(as.character(svm_pred),svm_prob) #need to chang svm_pred to numeric for agreescore creation
        colnames(results_table)[1] <- "SVM_LABEL"
        colnames(results_table)[2] <- "SVM_PROB"
    } else
    
    if (pmatch("slda",class(model),nomatch=0) > 0){
        slda_results <- predict(model,data.frame(as.matrix(container@classification_matrix)),...)
		slda_pred <- apply(slda_results$posterior,1,extract_label_from_prob_names) #Extract Label Based on Probability
		slda_prob <- apply(slda_results$posterior,1,extract_maximum_prob) #Extract Highest Probability
        
        results_table <- data.frame(as.character(slda_pred),slda_prob)
        colnames(results_table)[1] <- "SLDA_LABEL"
		colnames(results_table)[2] <- "SLDA_PROB"
    } else

    if (pmatch("LogitBoost",class(model),nomatch=0) > 0) {
        lboost_results <- predict(model,xtest=as.matrix(container@classification_matrix),type="raw",...) #Probability
        lboost_pred <- apply(lboost_results,1,extract_label_from_prob_names) #Extract Label Based on Probability
        lboost_prob <- apply(lboost_results,1,extract_maximum_prob) #Extract Highest Probability
        
        results_table <- data.frame(as.character(lboost_pred),lboost_prob)
        colnames(results_table)[1] <- "LOGITBOOST_LABEL"
        colnames(results_table)[2] <- "LOGITBOOST_PROB"
    } else
    
    if (pmatch("classbagg",class(model),nomatch=0) > 0) {
        bagging_results <- predict(model,newdata=data.frame(as.matrix(container@classification_matrix)), type=c("prob"),...)
        bagging_pred <- apply(bagging_results,1,extract_label_from_prob_names) #Extract Label Based on Probability
        bagging_prob <- apply(bagging_results,1,extract_maximum_prob) 
        
        results_table <- data.frame(as.character(bagging_pred),bagging_prob)
        colnames(results_table)[1] <- "BAGGING_LABEL"
        colnames(results_table)[2] <- "BAGGING_PROB"
    } else
    
    if (pmatch("randomForest",class(model),nomatch=0) > 0){
        rf_results <- predict(model,newdata=as.matrix(container@classification_matrix),type="prob",...)
		rf_pred <- apply(rf_results,1,extract_label_from_prob_names)
        rf_prob <- apply(rf_results,1,extract_maximum_prob)

        results_table <- data.frame(as.character(rf_pred),rf_prob)
        colnames(results_table)[1] <- "FORESTS_LABEL"
        colnames(results_table)[2] <- "FORESTS_PROB"
    } else
    
    if (pmatch("glmnet",class(model),nomatch=0) > 0){
		classification_matrix <- as(as.matrix.csc(container@classification_matrix),"dgCMatrix")
		#colnames(classification_matrix) <- container@column_names
        glmnet_results <- predict(model,newx=classification_matrix,s=s,type="response",...)
        glmnet_pred <- apply(glmnet_results[,,1],1,extract_label_from_prob_names) 
        glmnet_prob <- apply(glmnet_results,1,extract_maximum_prob) 
        
        results_table <- data.frame(as.character(glmnet_pred),glmnet_prob)
        colnames(results_table)[1] <- "GLMNET_LABEL"
        colnames(results_table)[2] <- "GLMNET_PROB"
    } else
    
    if (pmatch("tree",class(model),nomatch=0) > 0){
        tree_results <- predict(model,newdata=data.frame(as.matrix(container@classification_matrix)), type="vector",...)
        tree_pred <- apply(tree_results,1,extract_label_from_prob_names)
        tree_prob <- apply(tree_results,1,extract_maximum_prob) 
        
        results_table <- data.frame(as.character(tree_pred),tree_prob)
        colnames(results_table)[1] <- "TREE_LABEL"
        colnames(results_table)[2] <- "TREE_PROB"
    } else

    if (pmatch("nnet",class(model),nomatch=0) > 0){
        nnet_results <- predict(model,newdata=data.frame(as.matrix(container@classification_matrix)),...) #probabilities
        nnet_pred <- apply(nnet_results,1,extract_label_from_prob_names) #Extract Highest Probability Score
        nnet_prob <- apply(nnet_results,1,extract_maximum_prob) #Extract Probability
        
        results_table <- data.frame(as.character(nnet_pred),nnet_prob)
        colnames(results_table)[1] <- "NNETWORK_LABEL"
        colnames(results_table)[2] <- "NNETWORK_PROB"
    } else
							   
	if (pmatch("maxent",class(model),nomatch=0) > 0) {
		maxent_results <- predict(model,container@classification_matrix,...)
		maxent_pred <- maxent_results[,1]
		maxent_prob <- apply(maxent_results[,-1],1,extract_maximum_prob)
		
		results_table <- data.frame(as.character(maxent_pred),as.vector(maxent_prob,mode="numeric"))
		colnames(results_table)[1] <- "MAXENTROPY_LABEL"
		colnames(results_table)[2] <- "MAXENTROPY_PROB"
	}
	
	return(results_table)
}

