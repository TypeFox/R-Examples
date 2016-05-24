# LOAD THE PACKAGE
library(RTextTools)

# SET THE SEED AND LOAD THE DATA
set.seed(95616)
data(USCongress)

# CREATE THE DOCUMENT-TERM MATRIX AND WRAP THE DATA IN A CONTAINER
doc_matrix <- create_matrix(USCongress$text, language="english", removeNumbers=TRUE, stemWords=TRUE, removeSparseTerms=.998)
container <- create_container(doc_matrix, USCongress$major, trainSize=1:4000, testSize=4001:4449, virgin=FALSE)

# TRAIN THE ALGORITHMS USING THE CONTAINER
# ALTERNATIVELY, train_models(container, c("SVM","GLMNET","MAXENT","SLDA","BOOSTING","BAGGING","RF","NNET","TREE"))
SVM <- train_model(container,"SVM")
GLMNET <- train_model(container,"GLMNET")
MAXENT <- train_model(container,"MAXENT")
SLDA <- train_model(container,"SLDA")
BOOSTING <- train_model(container,"BOOSTING")
BAGGING <- train_model(container,"BAGGING")
RF <- train_model(container,"RF")
NNET <- train_model(container,"NNET")
TREE <- train_model(container,"TREE")

# CLASSIFY THE TESTING DATA USING THE TRAINED MODELS.
# ALTERNATIVELY, classify_models(container, list_of_trained_models)
SVM_CLASSIFY <- classify_model(container, SVM)
GLMNET_CLASSIFY <- classify_model(container, GLMNET)
MAXENT_CLASSIFY <- classify_model(container, MAXENT)
SLDA_CLASSIFY <- classify_model(container, SLDA)
BOOSTING_CLASSIFY <- classify_model(container, BOOSTING)
BAGGING_CLASSIFY <- classify_model(container, BAGGING)
RF_CLASSIFY <- classify_model(container, RF)
NNET_CLASSIFY <- classify_model(container, NNET)
TREE_CLASSIFY <- classify_model(container, TREE)

# CREATE THE ANALYTICS USING THE RESULTS FROM ALL THE ALGORITHMS
analytics <- create_analytics(container,cbind(SVM_CLASSIFY, SLDA_CLASSIFY, 
	BOOSTING_CLASSIFY, BAGGING_CLASSIFY, RF_CLASSIFY, GLMNET_CLASSIFY, 
	NNET_CLASSIFY, TREE_CLASSIFY, MAXENT_CLASSIFY))

# DEMONSTRATION OF HOW TO WRITE THE DATA OUT TO A .CSV FILE
# write.csv(analytics@document_summary,"DocumentSummary.csv")
# write.csv(analytics@topic_summary,"TopicSummary.csv")
# write.csv(analytics@algorithm_summary,"AlgorithmSummary.csv")
# write.csv(analytics@ensemble_summary,"EnsembleSummary.csv")

