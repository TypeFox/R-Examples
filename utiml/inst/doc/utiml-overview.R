## ------------------------------------------------------------------------
library("utiml")

## ---- echo=FALSE, results='asis'-----------------------------------------
bl <- data.frame(
  Use = c("CART", "C5.0", "J48", "KNN", "MAJORITY", "NB", "RANDOM", "RF", "SVM"),
  Name = c("Classification and regression trees", "C5.0 Decision Trees and Rule-Based Models", "Java implementation of the C4.5", "K Nearest Neighbor", "Majority class prediction", "Naive Bayes", "Random prediction", "Random Forest", "Support Vector Machine"),
  Package = c("rpart", "C50", "RWeka and rJava", "kknn", "-", "e1071", "-", "randomForest", "e1071"),
  Call = c("rpart::rpart(...)", "C50::C5.0(...)", "RWeka::J48(...)", "kknn::kknn(...)", "-", "e1071::naiveBayes(...)", "-", "randomForest::randomForest(...)", "e1071::svm(...)")
)
knitr::kable(bl)

## ---- echo=FALSE, results='asis'-----------------------------------------
mts <- data.frame(
  Method = c("br", "brplus", "cc", "ctrl", "dbr", "ebr", "ecc", "mbr", "ns", "prudent", "rdbr"),
  Name = c("Binary Relevance (BR)", "BR+", "Classifier Chains", "ConTRolled Label correlation exploitation (CTRL)",
           "Dependent Binary Relevance (DBR)", "Ensemble of Binary Relevance (EBR)", "Ensemble of Classifier Chains (ECC)", "Meta-Binary Relevance (MBR or 2BR)", "Nested Stacking (NS)", "Pruned and Confident Stacking Approach (Prudent)", "Recursive Dependent Binary Relevance (RDBR)"),
  Approach = c("one-agains-all", "one-agains-all; stacking", "one-agains-all; stacking", "one-agains-all; binary-ensemble", "one-agains-all; stacking", "one-agains-all; ensemble", "one-agains-all; ensemble; stacking",
"one-agains-all; stacking", "one-agains-all; stacking", "one-agains-all; binary-ensemble; stacking", "one-agains-all; stacking")
)
knitr::kable(mts)

## ------------------------------------------------------------------------
toy <- create_holdout_partition(toyml)
brmodel <- br(toy$train, "SVM")
prediction <- predict(brmodel, toy$test)

# Using the test dataset and the prediction
result <- multilabel_evaluate(toy$test, prediction)
print(round(result, 3))

# Build a confusion matrix
confmat <- multilabel_confusion_matrix(toy$test, prediction)
result <- multilabel_evaluate(confmat)
print(confmat)

## ------------------------------------------------------------------------
# Example-based measures
result <- multilabel_evaluate(confmat, "example-based")
print(names(result))

# Subset accuracy, F1 measure and hamming-loss
result <- multilabel_evaluate(confmat, c("subset-accuracy", "F1", "hamming-loss"))
print(names(result))

# Ranking and label-basedd measures
result <- multilabel_evaluate(confmat, c("label-based", "ranking"))
print(names(result))

# To see all the supported measures you can try
multilabel_measures()

## ---- echo=FALSE, results='asis'-----------------------------------------
## 8. How to extend utiml

### 8.1 Create a new Multi-label Method

### 8.2 Create a new base Learner

