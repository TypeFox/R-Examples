### R code from vignette source 'pensim.Rnw'

###################################################
### code chunk number 1: foo
###################################################
options(keep.source = TRUE)
options(width = 60)
options(eps = FALSE)
foo <- packageDescription("pensim")


###################################################
### code chunk number 2: dataprep
###################################################
library(pensim)
data(beer.exprs)
data(beer.survival)
##select just 100 genes to speed computation, just for the sake of example:
set.seed(1)
beer.exprs.sample <- beer.exprs[sample(1:nrow(beer.exprs), 100), ]
#
gene.quant <- apply(beer.exprs.sample, 1, quantile, probs=0.75)
dat.filt <- beer.exprs.sample[gene.quant>log2(100), ]
gene.iqr <- apply(dat.filt, 1, IQR)
dat.filt <- as.matrix(dat.filt[gene.iqr>0.5, ])
dat.filt <- t(dat.filt)
dat.filt <- data.frame(dat.filt)
#
library(survival)
surv.obj <- Surv(beer.survival$os, beer.survival$status)


###################################################
### code chunk number 3: lassotest
###################################################
library(penalized)
testfit <- optL1(response=surv.obj,
                 penalized=dat.filt,
                 fold=5,
                 maxlambda1=5,
                 positive=FALSE,
                 standardize=TRUE,
                 trace=FALSE)


###################################################
### code chunk number 4: opt.nested.crossval
###################################################
set.seed(1)
preds <- opt.nested.crossval(outerfold=5, nprocessors=1, #opt.nested.crossval arguments
             optFUN="opt1D", scaling=FALSE,                     #opt.splitval arguments
             setpen="L1", nsim=1,                                      #opt1D arguments
             response=surv.obj,                    #rest are penalized::optl1 arguments
             penalized=dat.filt,
             fold=5,
             positive=FALSE,
             standardize=TRUE,
             trace=FALSE)


###################################################
### code chunk number 5: coxfit
###################################################
coxfit.continuous <- coxph(surv.obj~preds)
summary(coxfit.continuous)


###################################################
### code chunk number 6: dichot
###################################################
preds.dichot <- preds > median(preds)


###################################################
### code chunk number 7: prelim
###################################################
options(device="pdf")


###################################################
### code chunk number 8: ROCplot
###################################################
if(require(survivalROC)){
  nobs <- length(preds)
  cutoff <- 12
  preds.roc <- survivalROC(Stime=beer.survival$os, status=beer.survival$status,
                         marker=preds, predict.time=cutoff, span = 0.01*nobs^(-0.20))
  plot(preds.roc$FP, preds.roc$TP, type="l", xlim=c(0, 1), ylim=c(0, 1),   
       xlab=paste( "FP", "\n", "AUC = ", round(preds.roc$AUC, 3)), 
       lty=2,
       ylab="TP", main="LASSO predictions\n ROC curve at 12 months")
  abline(0, 1)
}


###################################################
### code chunk number 9: full.model
###################################################
beer.coefs <- opt1D(setpen="L1",nsim=1,                        
             response=surv.obj,                         
             penalized=dat.filt,
             fold=5,
             maxlambda1=5,
             positive=FALSE,
             standardize=TRUE,
             trace=FALSE)


###################################################
### code chunk number 10: unpenalized.eg
###################################################
beer.coefs.unpen <- opt1D(setpen="L1", nsim=1,                        
             response=surv.obj,                         
             penalized=dat.filt[-1],   # This is equivalent to dat.filt[,-1]
             unpenalized=dat.filt[1],  
             fold=5,
             maxlambda1=5,
             positive=FALSE,
             standardize=TRUE,
             trace=FALSE)


###################################################
### code chunk number 11: lookatcoefs
###################################################
beer.coefs[1, 1:5]        #example output with no unpenalized covariates
beer.coefs.unpen[1, 1:5]  #example output with first covariate unpenalized


###################################################
### code chunk number 12: genbinary
###################################################
set.seed(9)
x <- create.data(nvars=c(15, 5),
                 cors=c(0, 0.8),
                 associations=c(0, 2),
                 firstonly=c(TRUE, TRUE),
                 nsamples=50,
                 response="binary",
                 logisticintercept=0.5)


###################################################
### code chunk number 13: lookbinary
###################################################
summary(x)
x$summary


###################################################
### code chunk number 14: fitmodel
###################################################
simplemodel <- glm(outcome ~ ., data=x$data, family=binomial)
summary(simplemodel)


###################################################
### code chunk number 15: binarylassodemo
###################################################
lassofit <- opt1D(nsim=3, nprocessors=1, setpen="L1", penalized=x$data[1:20],
                 response=x$data[, "outcome"], trace=FALSE, fold=10)
print(lassofit)


###################################################
### code chunk number 16: heatmap
###################################################
dat <- t(as.matrix(x$data[, -match("outcome", colnames(x$data))]))
heatmap(dat, ColSideColors=ifelse(x$data$outcome==0, "black", "white"))


###################################################
### code chunk number 17: survoutcome
###################################################
set.seed(1)
x <- create.data(nvars=c(15, 5),
                 cors=c(0, 0.8),
                 associations=c(0, 0.5),
                 firstonly=c(TRUE, TRUE),
                 nsamples=50,
                 censoring=c(2, 10),
                 response="timetoevent")


###################################################
### code chunk number 18: howmanycensored
###################################################
sum(x$data$cens==0)/nrow(x$data)


###################################################
### code chunk number 19: simulatedKM
###################################################
library(survival)
surv.obj <- Surv(x$data$time, x$data$cens)
plot(survfit(surv.obj~1), ylab="Survival probability", xlab="time")


###################################################
### code chunk number 20: pensim.Rnw:365-366
###################################################
sessionInfo()


