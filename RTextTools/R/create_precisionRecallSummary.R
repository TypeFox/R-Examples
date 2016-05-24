create_precisionRecallSummary <- function(container, classification_results, b_value=1) {	
	confusion <- function(true,pred) {
		conf_out <- table(factor(true,levels=sort(unique(true))),factor(pred,levels=sort(unique(true))))
		return(conf_out)
	}
	
	precision <- function(confusion) {
		correct <- diag(confusion)
		precision_sums <- colSums(confusion)
		precision_out <- round (correct/precision_sums,2)
		return(precision_out)
	}
	
	recall <- function(confusion) {
		correct <- diag(confusion)
		recall_sums <- rowSums(confusion)
		recall_out <- round (correct/recall_sums,2)
		return(recall_out)
	}
	
# takes v_value set to 1 weighs precision/recall equally; vector of precision values, vector recall values
	fscore <- function(b_value,precision,recall){
		B <- b_value
		fscore <- ((B^2+1) * precision * recall) / ((B^2 * precision) + recall)
		return(fscore)
	}
	
	fscores_out <- function(b_value,precision,recall) {
		fscores_out <- NULL
		for (i in seq_along(precision)) {
			fscores_out[i] <- round(fscore(b_value,precision[i],recall[i]),2)
		}
		return(fscores_out)
	}
	
	scores <- create_scoreSummary(container, classification_results)
	
	true <- container@testing_codes
	columns <- colnames(scores)
	results <- c()
	
	if (pmatch("SVM_LABEL",columns,nomatch=0) > 0) {
		pred <- scores$SVM_LABEL
		conf <- confusion(true,pred)
		
		svm_precision <- precision(conf)
		
		svm_recall <- recall(conf)
		
		svm_fscore <- fscores_out (b_value,svm_precision,svm_recall)
		svm_results <- cbind(SVM_PRECISION=svm_precision,SVM_RECALL=svm_recall,SVM_FSCORE=svm_fscore)
		
		results <- cbind(results,svm_results)
	}
	
	if (pmatch("SLDA_LABEL",columns,nomatch=0) > 0) {
		pred <- scores$SLDA_LABEL
		conf <- confusion(true,pred)
		
		slda_precision <- precision(conf)
		
		slda_recall <- recall(conf)
		
		slda_fscore <- fscores_out (b_value,slda_precision,slda_recall)
		slda_results <- cbind(SLDA_PRECISION=slda_precision,SLDA_RECALL=slda_recall,SLDA_FSCORE=slda_fscore)
		
		results <- cbind(results,slda_results)
	}
	
	if (pmatch("LOGITBOOST_LABEL",columns,nomatch=0) > 0) {
		pred <- scores$LOGITBOOST_LABEL
		conf <- confusion(true,pred)
		
		boosting_precision <- precision(conf)
		
		boosting_recall <- recall(conf)
		
		boosting_fscore <- fscores_out (b_value,boosting_precision,boosting_recall)
		boosting_results <- cbind(LOGITBOOST_PRECISION=boosting_precision,LOGITBOOST_RECALL=boosting_recall,LOGITBOOST_FSCORE=boosting_fscore)
		
		results <- cbind(results,boosting_results)
	}
	
	if (pmatch("BAGGING_LABEL",columns,nomatch=0) > 0) {
		pred <- scores$BAGGING_LABEL
		conf <- confusion(true,pred)
		
		bagging_precision <- precision(conf)
		
		bagging_recall <- recall(conf)
		
		bagging_fscore <- fscores_out (b_value,bagging_precision,bagging_recall)
		bagging_results <- cbind(BAGGING_PRECISION=bagging_precision,BAGGING_RECALL=bagging_recall,BAGGING_FSCORE=bagging_fscore)
		
		results <- cbind(results,bagging_results)
	}
	
	if (pmatch("FORESTS_LABEL",columns,nomatch=0) > 0) {
		pred <- scores$FORESTS_LABEL
		conf <- confusion(true,pred)
		
		rf_precision <- precision(conf)
		
		rf_recall <- recall(conf)
		
		rf_fscore <- fscores_out (b_value,rf_precision,rf_recall)
		rf_results <- cbind(FORESTS_PRECISION=rf_precision,FORESTS_RECALL=rf_recall,FORESTS_FSCORE=rf_fscore)
		
		results <- cbind(results,rf_results)
	}
	
	if (pmatch("GLMNET_LABEL",columns,nomatch=0) > 0) {
		pred <- scores$GLMNET_LABEL
		conf <- confusion(true,pred)
		
		glmnet_precision <- precision(conf)
		
		glmnet_recall <- recall(conf)
		
		glmnet_fscore <- fscores_out (b_value,glmnet_precision,glmnet_recall)
		glmnet_results <- cbind(GLMNET_PRECISION=glmnet_precision,GLMNET_RECALL=glmnet_recall,GLMNET_FSCORE=glmnet_fscore)
		
		results <- cbind(results,glmnet_results)
	}
	
	if (pmatch("TREE_LABEL",columns,nomatch=0) > 0) {
		pred <- scores$TREE_LABEL
		conf <- confusion(true,pred)
		
		tree_precision <- precision(conf)
		
		tree_recall <- recall(conf)
		
		tree_fscore <- fscores_out (b_value,tree_precision,tree_recall)
		tree_results <- cbind(TREE_PRECISION=tree_precision,TREE_RECALL=tree_recall,TREE_FSCORE=tree_fscore)
		
		results <- cbind(results,tree_results)
	}
	
	if (pmatch("NNETWORK_LABEL",columns,nomatch=0) > 0) {
		pred <- scores$NNETWORK_LABEL
		conf <- confusion(true,pred)
		
		nnet_precision <- precision(conf)
		
		nnet_recall <- recall(conf)
		
		nnet_fscore <- fscores_out (b_value,nnet_precision,nnet_recall)
		nnet_results <- cbind(NNETWORK_PRECISION=nnet_precision,NNETWORK_RECALL=nnet_recall,NNETWORK_FSCORE=nnet_fscore)
		
		results <- cbind(results,nnet_results)
	}
	
	if (pmatch("MAXENTROPY_LABEL",columns,nomatch=0) > 0) {
		pred <- scores$MAXENTROPY_LABEL
		conf <- confusion(true,pred)
		
		maxent_precision <- precision(conf)
		
		maxent_recall <- recall(conf)
		
		maxent_fscore <- fscores_out (b_value,maxent_precision,maxent_recall)
		maxent_results <- cbind(MAXENTROPY_PRECISION=maxent_precision,MAXENTROPY_RECALL=maxent_recall,MAXENTROPY_FSCORE=maxent_fscore)
		
		results <- cbind(results,maxent_results)
	}
	
	return(results)
}