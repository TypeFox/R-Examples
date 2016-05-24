library(partDSA)
library(MASS)

set.seed(6442)

n <- nrow(Boston)
tr.n <- floor(n / 2)
ts.n <- n - tr.n

train.index <- sample(1:n, tr.n, replace=FALSE)
test.index <- c(1:n)[-train.index]
vars <- Boston[,-14]
outcome <- Boston[,14]

x <- vars[train.index,]
y <- outcome[train.index]
wt <- rep(1, tr.n)

x.test <- vars[test.index,]
y.test <- outcome[test.index]
wt.test <- rep(1, ts.n)

control <- DSA.control(vfold=1)  # disable cross-validation
results <- partDSA(x, y, wt, x.test, y.test, wt.test, control=control)
print(results)

# View the "model" object using the Java visualization program
showDSA(results)
