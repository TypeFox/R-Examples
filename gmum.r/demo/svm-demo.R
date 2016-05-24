library(gmum.r)
library(caret)

# Load a dataset, here we have provided an example 
data(svm_breast_cancer_dataset)
ds <- svm.breastcancer.dataset

# Create CV folds
K <- 5

folds <- createFolds(ds$X1, k=K)
mean_acc <- 0

# SVM model needs to know how the labels depend on data
formula <- X1~. 

# Iterate through folds
for ( i in seq(1,K,1) ) {
  
  # Get training and testing data 
  train <- ds[-folds[[i]],]
  test <- ds[folds[[i]],]
  
  # Train SVM model
  svm <- SVM(formula, train, lib="libsvm", kernel="linear", prep = "2e", C=10);
  
  # Plot one of the SVMs using PCA
  if (i == 1) plot(svm, mode="pca")
  
  # Seperate lables in test data
  test_x <- subset(test, select = -c(X1))
  target <- test[,"X1"]
  
  # predict on test data
  pred <- predict(svm, test_x)
  
  # calculate classification accuracy
  acc <- svm.accuracy(prediction=pred, target=target)
  mean_acc <- mean_acc + acc  
}

# Display mean accuracy
print(sprintf("mean SVM accuracy after %i folds: %f ", K, mean_acc/K))

# Print short summray of the last trained svm
summary(svm)
