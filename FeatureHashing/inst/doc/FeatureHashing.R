## ----setup, include=FALSE------------------------------------------------
# library(knitcitations)
# bib <- read.bibtex("README.bib")
# citep(bib[[1]])

## ----lr------------------------------------------------------------------
# The following script assumes that the data.frame
# of the training dataset and testing dataset are 
# assigned to variable `ipinyou.train` and `ipinyou.test`
# respectively

library(FeatureHashing)

# Checking version.
stopifnot(packageVersion("FeatureHashing") >= package_version("0.9"))

data(ipinyou)
f <- ~ IP + Region + City + AdExchange + Domain +
  URL + AdSlotId + AdSlotWidth + AdSlotHeight +
  AdSlotVisibility + AdSlotFormat + CreativeID +
  Adid + split(UserTag, delim = ",")
# if the version of FeatureHashing is 0.8, please use the following command:
# m.train <- as(hashed.model.matrix(f, ipinyou.train, 2^16, transpose = FALSE), "dgCMatrix")
m.train <- hashed.model.matrix(f, ipinyou.train, 2^16)
m.test <- hashed.model.matrix(f, ipinyou.test, 2^16)

# logistic regression with glmnet

library(glmnet)

cv.g.lr <- cv.glmnet(m.train, ipinyou.train$IsClick,
  family = "binomial")#, type.measure = "auc")
p.lr <- predict(cv.g.lr, m.test, s="lambda.min")
auc(ipinyou.test$IsClick, p.lr)

## ----xgboost-------------------------------------------------------------
# GBDT with xgboost
if(require("xgboost")){
  cv.g.gdbt <- xgboost(m.train, ipinyou.train$IsClick, max.depth=7, eta=0.1, subsample = 0.7, colsample_bytree = 0.7,
    nround = 100, objective = "binary:logistic", verbose = ifelse(interactive(), 1, 0))
  p.lm <- predict(cv.g.gdbt, m.test)
  glmnet::auc(ipinyou.test$IsClick, p.lm)  
}


## ----ftprl---------------------------------------------------------------
source(system.file("ftprl.R", package = "FeatureHashing"))

m.train <- hashed.model.matrix(f, ipinyou.train, 2^16, transpose = TRUE)
ftprl <- initialize.ftprl(0.1, 1, 0.1, 0.1, 2^16)
ftprl <- update.ftprl(ftprl, m.train, ipinyou.train$IsClick, predict = TRUE)
auc(ipinyou.train$IsClick, attr(ftprl, "predict"))

