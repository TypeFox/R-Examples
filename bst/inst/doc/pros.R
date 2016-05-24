### R code from vignette source 'pros.Rnw'

###################################################
### code chunk number 1: pros.Rnw:35-36
###################################################
options(prompt = "R> ", continue = " ", width = 70, digits =4, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: pros.Rnw:38-50
###################################################
library("bst")
library("pROC")
### benign prostate hyperplasia
bph <- read.table(file=system.file("extdata", "BPH.txt", package="bst"))
### early stage prostate cancer, 
ccd <- read.table(file=system.file("extdata", "CCD.txt", package="bst"))
### late stage prostate cancer
cab <- read.table(file=system.file("extdata", "CAB.txt", package="bst"))
### age-matched healthy men
hm <- read.table(file=system.file("extdata", "control.txt", package="bst"))
### mass spectrometry protein
mz <- read.table(file=system.file("extdata", "prostate_mz.txt", package="bst"))


###################################################
### code chunk number 3: pros.Rnw:54-56
###################################################
dat <- t(cbind(hm, bph, cab, ccd))
y <- c(rep(-1, dim(hm)[2] + dim(bph)[2]), rep(1,  dim(cab)[2] + dim(ccd)[2]))


###################################################
### code chunk number 4: pros.Rnw:60-69
###################################################
myuniq <- function(x){
  if(sum(x == 0) <= length(x) * 0.95)
  return(TRUE)
  else return(FALSE)
}
res <- apply(dat, 2, myuniq)
dat <- dat[,res]
rownames(dat) <- NULL
colnames(dat) <- mz[res,1]


###################################################
### code chunk number 5: pros.Rnw:73-79
###################################################
ntrain <- floor(length(y)*0.75)
set.seed(13)
q <- sample(length(y))
y.tr <- y[q][1:ntrain]; y.te <- y[q][-(1:ntrain)]
X <- dat[q,]
x.tr <- X[1:ntrain,];  x.te <- X[-(1:ntrain),] 


###################################################
### code chunk number 6: pros.Rnw:82-90
###################################################
dat.m1 <- bst(x=x.tr, y=y.tr, ctrl = bst_control(mstop=400), family = "hinge2")
pred <- predict(dat.m1, x.te)
### misclassification error
mean(abs(y.te-sign(pred))/2)
### area under the curve (AUC) of receiver operating characteristic (ROC)
auc(y.te, pred)
### number of variables selected
length(dat.m1$xselect)


###################################################
### code chunk number 7: pros.Rnw:93-103
###################################################
dat.m2 <- bst(x=x.tr, y=y.tr, family="hinge2", ctrl = bst_control(mstop=500, 
twinboost=TRUE, twintype=2, coefir=coef(dat.m1), f.init=predict(dat.m1), 
xselect.init = dat.m1$xselect))
pred <- predict(dat.m2, x.te)
### misclassification error
mean(abs(y.te-sign(pred))/2)
### AUC with twin boosting
auc(y.te, pred)
### number of variables selected
length(dat.m2$xselect)


###################################################
### code chunk number 8: pros.Rnw:106-118
###################################################
dat <- t(cbind(cab, ccd, hm))
y <- c(rep(1, dim(cab)[2] + dim(ccd)[2]), rep(-1, dim(hm)[2]))
res <- apply(dat, 2, myuniq)
dat <- dat[,res]
rownames(dat) <- NULL
colnames(dat) <- mz[res,1]
ntrain <- floor(length(y)*0.75)
set.seed(13)
q <- sample(length(y))
y.tr <- y[q][1:ntrain]; y.te <- y[q][-(1:ntrain)]
X <- dat[q,]
x.tr <- X[1:ntrain,];  x.te <- X[-(1:ntrain),] 


###################################################
### code chunk number 9: pros.Rnw:121-129
###################################################
dat.m1 <- bst(x=x.tr, y=y.tr, ctrl = bst_control(mstop=400), family = "hinge2")
pred <- predict(dat.m1, x.te)
### misclassification error
mean(abs(y.te-sign(pred))/2)
### AUC
auc(y.te, pred)
### number of variables selected
length(dat.m1$xselect)


###################################################
### code chunk number 10: pros.Rnw:132-142
###################################################
dat.m2 <- bst(x=x.tr, y=y.tr, family="hinge2", ctrl = bst_control(mstop=200, 
twinboost=TRUE, twintype=2, coefir=coef(dat.m1), f.init=predict(dat.m1), 
xselect.init = dat.m1$xselect))
pred <- predict(dat.m2, x.te)
### misclassification error
mean(abs(y.te-sign(pred))/2)
### AUC with twin boosting
auc(y.te, pred)
### number of variables selected
length(dat.m2$xselect)


###################################################
### code chunk number 11: pros.Rnw:145-157
###################################################
dat <- t(cbind(bph, cab, ccd))
y <- c(rep(-1, dim(bph)[2]), rep(1,  dim(cab)[2] + dim(ccd)[2]))
res <- apply(dat, 2, myuniq)
dat <- dat[,res]
rownames(dat) <- NULL
colnames(dat) <- mz[res,1]
ntrain <- floor(length(y)*0.75)
set.seed(13)
q <- sample(length(y))
y.tr <- y[q][1:ntrain]; y.te <- y[q][-(1:ntrain)]
X <- dat[q,]
x.tr <- X[1:ntrain,];  x.te <- X[-(1:ntrain),] 


###################################################
### code chunk number 12: pros.Rnw:160-168
###################################################
dat.m1 <- bst(x=x.tr, y=y.tr, ctrl = bst_control(mstop=400), family = "hinge2")
pred <- predict(dat.m1, x.te)
### misclassification error
mean(abs(y.te-sign(pred))/2)
### AUC
auc(y.te, pred)
### number of variables selected
length(dat.m1$xselect)


###################################################
### code chunk number 13: pros.Rnw:171-181
###################################################
dat.m2 <- bst(x=x.tr, y=y.tr, family="hinge2", ctrl = bst_control(mstop=500, 
twinboost=TRUE, twintype=2, coefir=coef(dat.m1), f.init=predict(dat.m1), 
xselect.init = dat.m1$xselect))
pred <- predict(dat.m2, x.te)
### misclassification error
mean(abs(y.te-sign(pred))/2)
### AUC with twin boosting
auc(y.te, pred)
### number of variables selected
length(dat.m2$xselect)


###################################################
### code chunk number 14: pros.Rnw:184-196
###################################################
dat <- t(cbind(bph, hm))
y <- c(rep(1, dim(bph)[2]), rep(-1, dim(hm)[2]))
res <- apply(dat, 2, myuniq)
dat <- dat[,res]
rownames(dat) <- NULL
colnames(dat) <- mz[res,1]
ntrain <- floor(length(y)*0.75)
set.seed(13)
q <- sample(length(y))
y.tr <- y[q][1:ntrain]; y.te <- y[q][-(1:ntrain)]
X <- dat[q,]
x.tr <- X[1:ntrain,];  x.te <- X[-(1:ntrain),] 


###################################################
### code chunk number 15: pros.Rnw:199-207
###################################################
dat.m1 <- bst(x=x.tr, y=y.tr, ctrl = bst_control(mstop=400), family = "hinge2")
pred <- predict(dat.m1, x.te)
### misclassification error
mean(abs(y.te-sign(pred))/2)
### AUC
auc(y.te, pred)
### number of variables selected
length(dat.m1$xselect)


###################################################
### code chunk number 16: pros.Rnw:210-220
###################################################
dat.m2 <- bst(x=x.tr, y=y.tr, family="hinge2", ctrl = bst_control(mstop=500, 
twinboost=TRUE, twintype=2, coefir=coef(dat.m1), f.init=predict(dat.m1), 
xselect.init = dat.m1$xselect))
pred <- predict(dat.m2, x.te)
### misclassification error
mean(abs(y.te-sign(pred))/2)
### AUC with twin boosting
auc(y.te, pred)
### number of variables selected
length(dat.m2$xselect)


