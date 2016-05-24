### R code from vignette source 'weibull.Rnw'

###################################################
### code chunk number 1: init
###################################################
if(require(SurvRegCensCov) == FALSE){install.packages("SurvRegCensCov"); library("SurvRegCensCov")}


###################################################
### code chunk number 2: source
###################################################
library(survival)
data(larynx)


###################################################
### code chunk number 3: WeibullReg
###################################################
WeibullReg(Surv(time, death) ~ factor(stage) + age, data=larynx)


###################################################
### code chunk number 4: coxph
###################################################
cph <- coxph(Surv(time, death) ~ factor(stage) + age, data=larynx)
summary(cph)$conf.int


###################################################
### code chunk number 5: DiagPlot
###################################################
WeibullDiag(Surv(time, death) ~ factor(stage), data = larynx, 
            labels=c("Stage I", "Stage II", "Stage III", "Stage IV"))


