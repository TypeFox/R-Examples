# LOAD THE PACKAGE
library(RTextTools)

# SET THE SEED AND LOAD THE DATA
set.seed(95616)
data(USCongress)

# CREATE THE DOCUMENT-TERM MATRIX AND WRAP THE DATA IN A CONTAINER
doc_matrix <- create_matrix(USCongress$text, language="english", removeNumbers=TRUE, stemWords=TRUE, removeSparseTerms=.998)
container <- create_container(doc_matrix, USCongress$major, trainSize=1:4000, testSize=4001:4449, virgin=FALSE)

# DEMONSTRATION OF CROSS-VALIDATION
SVM <- cross_validate(container,4,"SVM")
MAXENT <- cross_validate(container,4,"MAXENT")
GLMNET <- cross_validate(container,4,"GLMNET")
SLDA <- cross_validate(container,4,"SLDA")
BAGGING <- cross_validate(container,4,"BAGGING")
BOOSTING <- cross_validate(container,4,"BOOSTING")
RF <- cross_validate(container,4,"RF")
NNET <- cross_validate(container,4,"NNET")
TREE <- cross_validate(container,4,"TREE")