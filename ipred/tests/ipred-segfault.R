library("ipred")
library("mlbench")
library("MASS")
library("survival")

actversion <- paste(R.version$major, R.version$minor, sep=".")
thisversion <- "1.7.0"

#if (compareVersion(actversion, thisversion) >= 0) {
#  RNGversion("1.6.2")
#}
set.seed(29081975)


# Classification

learn <- as.data.frame(mlbench.twonorm(200))
test <- as.data.frame(mlbench.twonorm(100))

# bagging

mod <- bagging(classes ~ ., data=learn, coob=TRUE, nbagg=10)
mod
predict(mod)[1:10]

# Double-Bagging

comb.lda <- list(list(model=lda, predict=function(obj, newdata)
                      predict(obj, newdata)$x))

mod <- bagging(classes ~ ., data=learn, comb=comb.lda, nbagg=10)
mod
predict(mod, newdata=test[1:10,])
predict(mod, newdata=test[1:10,], agg="aver")
predict(mod, newdata=test[1:10,], agg="wei")
predict(mod, newdata=test[1:10,], type="prob")
predict(mod, newdata=test[1:10,], type="prob", agg="aver")
predict(mod, newdata=test[1:10,], type="prob", agg="wei")

mypredict.lda <- function(object, newdata)
       predict(object, newdata = newdata)$class

errorest(classes ~ ., data=learn, model=lda, predict=mypredict.lda)
errorest(classes ~ ., data=learn, model=lda, predict=mypredict.lda,
  est.para=control.errorest(k=5, random=FALSE))

lapply(errorest(classes ~ ., data=learn, model=lda, predict=mypredict.lda,
  est.para=control.errorest(k=5, random=FALSE, getmodels=TRUE))$models, class)
errorest(classes ~ ., data=learn, model=bagging,
         est.para=control.errorest(k=2), nbagg=10)
errorest(classes ~ ., data=learn, model=bagging,
         est.para=control.errorest(k=2), nbagg=10, comb=comb.lda)
errorest(classes ~ ., data=learn, model=lda,
predict=mypredict.lda, estimator="boot")
errorest(classes ~ ., data=learn, model=lda,
predict=mypredict.lda, estimator="632plus")

# Regression

learn <- as.data.frame(mlbench.friedman1(100))
test <- as.data.frame(mlbench.friedman1(100))

# bagging

mod <- bagging(y ~ ., data=learn, coob=TRUE, nbagg=10)
mod
predict(mod)[1:10]

predict(mod, newdata=test[1:10,])
predict(mod, newdata=test[1:10,], agg="aver") 
predict(mod, newdata=test[1:10,], agg="wei")  
errorest(y ~ ., data=learn, model=lm)
errorest(y ~ ., data=learn, model=lm,
         est.para=control.errorest(k=5, random=FALSE))
lapply(errorest(y ~ ., data=learn, model=lm,
                est.para=control.errorest(k=5, random=FALSE, getmodels=TRUE))$models, class)
errorest(y ~ ., data=learn, model=lm, estimator="boot")

# survival

learn <- rsurv(100, model="C")
test <- rsurv(100, model="C")

mod <- bagging(Surv(time, cens) ~ ., data=learn, nbagg=10)
mod
predict(mod, newdata=test[1:10,])

#errorest(Surv(time, cens) ~ ., data=learn, model=bagging, 
#         est.para=list(k=2, random=FALSE), nbagg=5)
#errorest(Surv(time, cens) ~ ., data=learn, model=bagging, 
#         estimator="boot", nbagg=5, est.para=list(nboot=5))
#insert control.errorest
errorest(Surv(time, cens) ~ ., data=learn, model=bagging, 
         est.para=control.errorest(k=2, random=FALSE), nbagg=5)
errorest(Surv(time, cens) ~ ., data=learn, model=bagging, 
         estimator="boot", nbagg=5, est.para=control.errorest(nboot=5))

#lapply(errorest(Surv(time, cens) ~ ., data=learn, model=bagging, 
#         estimator="cv", nbagg=1, est.para=list(k=2, random=FALSE,
#         getmodels=TRUE))$models, class)
#insert control.errorest
lapply(errorest(Surv(time, cens) ~ ., data=learn, model=bagging, 
         estimator="cv", nbagg=1, est.para=control.errorest(k=2, random=FALSE,
         getmodels=TRUE))$models, class)

# bundling for regression

learn <- as.data.frame(mlbench.friedman1(100))
test <- as.data.frame(mlbench.friedman1(100))

comb <- list(list(model=lm, predict=predict.lm))

modc <- bagging(y ~ ., data=learn, nbagg=10, comb=comb)
modc
predict(modc, newdata=learn)[1:10]

# bundling for survival

while(FALSE) {
data("GBSG2", package = "ipred")
rcomb <- list(list(model=coxph, predict=predict))

mods <- bagging(Surv(time,cens) ~ ., data=GBSG2, nbagg=10, 
                comb=rcomb,  control=rpart.control(xval=0))
predict(mods, newdata=GBSG2[1:3,])

# test for method dispatch on integer valued responses
y <- sample(1:100, 100)
class(y)
x <- matrix(rnorm(100*5), ncol=5)
mydata <- as.data.frame(cbind(y, x))

cv(y, y ~ ., data=mydata, model=lm, predict=predict)
bootest(y, y ~ ., data=mydata, model=lm, predict=predict)
bagging(y ~., data=mydata, nbagg=10)
}
