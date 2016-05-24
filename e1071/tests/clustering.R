## cmeans clustering should also work on data frames

library(e1071)
data(iris)
set.seed(123)
cm1 <- cmeans(iris[,1:4], 10)
bc1 <- bclust(iris[,1:4], 3, base.centers=20,iter.base=50,
              base.method="cmeans")
