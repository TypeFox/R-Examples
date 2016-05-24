## Loading the data
library(supclust)
data(leukemia, package="supclust")

## A subsample that gave seg.faults earlier:
ii <- c(4:12, 14, 16:18, 20:22, 25, 27:31, 34, 36:38)

## Generating random test data
set.seed(724)
xN <- matrix(rnorm(750), 3, 250)

## Running time
.proctime00 <- proc.time()

## Fitting Wilma
fit1  <- wilma(leukemia.x, leukemia.y, noc = 1, trace = 2)
fit2  <- wilma(leukemia.x, leukemia.y, noc = 3, trace = 1)
fit3  <- wilma(leukemia.x, leukemia.y, noc = 4, once.per.clust = TRUE)
## quite *different* results on 32-bit and 64-bit Linux (Fedora 15, 2011):
## (why?)
fit4  <- wilma(leukemia.x, leukemia.y, noc = 5, flip= FALSE, trace = 1)

## Running time
cat('Time elapsed: ',proc.time() - .proctime00,'\n')

## Checking the output of fit1
fit1
summary(fit1)
plot(fit1)
fitted(fit1)

identical(predict(fit1), fitted(fit1))
predict(fit1, type = "cla")
predict(fit1, type = "fitt")

predict(fit1, newdata = xN)
predict(fit1, newdata = xN, type = "cla", classifier = "nnr")
predict(fit1, newdata = xN, type = "cla", classifier = "dlda")
predict(fit1, newdata = xN, type = "cla", classifier = "logreg")
predict(fit1, newdata = xN, type = "cla", classifier = "aggtrees")

## Checking the output of fit2
fit2
summary(fit2)
plot(fit2)
fitted(fit2)

identical(predict(fit2), fitted(fit2))
predict(fit2, type = "cla")
predict(fit2, type = "fitt")

predict(fit2, newdata = xN)
predict(fit2, newdata = xN, type = "cla", classifier = "nnr", noc = c(1,2,3))
predict(fit2, newdata = xN, type = "cla", classifier = "dlda", noc = c(1,3))
predict(fit2, newdata = xN, type = "cla", classifier = "logreg")
predict(fit2, newdata = xN, type = "cla", classifier = "aggtrees")

## Checking the output of fit3
fit3
summary(fit3)
plot(fit3)
fitted(fit3)

identical(predict(fit3), fitted(fit3))
predict(fit3, type = "cla")
predict(fit3, type = "fitt")

predict(fit3, newdata = xN)
predict(fit3, newdata = xN, type = "cla", classifier = "nnr", noc = c(1,2,3))
predict(fit3, newdata = xN, type = "cla", classifier = "dlda", noc = c(1,3))
predict(fit3, newdata = xN, type = "cla", classifier = "logreg")
predict(fit3, newdata = xN, type = "cla", classifier = "aggtrees")

## Checking the output of fit4  (different on differen platforms!)
fit4
summary(fit4)
plot(fit4)
fitted(fit4)

identical(predict(fit4), fitted(fit4))
predict(fit4, type = "cla")
predict(fit4, type = "fitt")

predict(fit4, newdata = xN)
predict(fit4, newdata = xN, type = "cla", classifier = "nnr", noc = c(1,2,3))
predict(fit4, newdata = xN, type = "cla", classifier = "dlda", noc = c(1,3))
predict(fit4, newdata = xN, type = "cla", classifier = "logreg")
predict(fit4, newdata = xN, type = "cla", classifier = "aggtrees")


## Running time
cat('Time elapsed: (total) ', proc.time() - .proctime00,'\n')
