library(fastAdaboost)
library(datasets)

#in case of iris, the first tree gets every sample right, and 
#hence we stop at the first tree
# should have adaboost error of 0.
data(iris)
set.seed(999)
iris_2 <- iris[iris$Species %in% c("setosa","versicolor"),]
iris_2$Species <- factor(iris_2$Species)
train_index <- sample(size=50,x=nrow(iris_2))
iris_2_train <- iris_2[train_index,]
iris_2_test <- iris_2[-train_index,]

test_that("Iris data works, for adaboost M1",{
  ada_obj <- adaboost(Species~., iris_2_train, 10)
  pred <- predict(ada_obj, iris_2_test)
  print(paste("Adaboost Error on iris:", pred$error))
})

test_that("Iris data works, for real adaboost ",{  
  ada_obj <- real_adaboost(Species~., iris_2_train, 10)
  pred <- predict(ada_obj, iris_2_test)
  print(paste("Real Adaboost Error on iris:", pred$error))
})
