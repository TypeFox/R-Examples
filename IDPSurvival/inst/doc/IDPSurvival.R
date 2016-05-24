### R code from vignette source 'IDPSurvival.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: idpsinst (eval = FALSE)
###################################################
## install.packages("IDPSurvival_VERSION.tar.gz",
##                  repos=NULL,type="source") ## from local file
## install.packages("IDPSurvival")  ## or from CRAN


###################################################
### code chunk number 2: <idps
###################################################
library("IDPSurvival")


###################################################
### code chunk number 3: IDPSurvival.Rnw:137-138
###################################################
options(width=60)


###################################################
### code chunk number 4: example
###################################################
n <- 30
lambda <- 5
X <- rexp(n, rate = lambda) # sample lifetimes
Y <- rexp(n, rate = lambda) # sample censoring times
status <- (X<Y)*1 
time <- X*status+Y*(1-status) 
dataset <- cbind(time,status)
dataset


###################################################
### code chunk number 5: isurvfit
###################################################
formula <- Surv(dataset[,1],dataset[,2]) ~ 1
fit <- isurvfit(formula, s=0.5, 
                conf.int=0.95,display=FALSE)
fit


###################################################
### code chunk number 6: isurvfit2
###################################################
dataset <- data.frame(time,status)
formula <- Surv(time,status) ~ 1
fit <- isurvfit(formula,dataset,s=0.5,display=FALSE)


###################################################
### code chunk number 7: mlekm
###################################################
plot(fit)
# Kaplan-Meier estimation
library(survival)
km <- survfit(formula,dataset)
lines(km,col='red')
legend('bottomleft',c("IDP","Kaplan-Meier"),lty=c(1,1),
       col=c('black','red'),pch=c('o','.'))


###################################################
### code chunk number 8: 2cov
###################################################
# Running isurvfit on lung (from survival package) with 
# two groups: Male and Female
data(lung,package='survival')
formula <- Surv(time,status) ~ sex
fit <- isurvfit(formula, lung)
legend('topright',c("Male","Female"),
       lty=c(1,1),col=c(1,2),pch=c(1,2))


###################################################
### code chunk number 9: 3cov
###################################################
# three groups: ph.ecog = 0, 1, 2
formula <- Surv(time,status) ~ ph.ecog
sel =!is.na(match(lung$ph.ecog,c(0,1,2)))
fit <- isurvfit(formula, lung, subset=sel)
legend('topright',names(fit$strata), 
       lty=rep(1,3),col=c(1:3),
       pch=c(1:3),title='ECOG performance score')


###################################################
### code chunk number 10: test
###################################################
# Tests for the lung cancer dataset if male are  
# more likely to live less than females
formula <- Surv(time,status) ~ sex
test <- isurvdiff(formula, lung, 
                  alternative='greater',
                  nsamples=100000)
print(test)


###################################################
### code chunk number 11: test2
###################################################
# Tests for the lung cancer dataset if male are more likely 
# to live less than females
formula <- Surv(time,status) ~ ph.ecog
test <- isurvdiff(formula, lung, groups=c(0,1),
                  alternative='two.sided', 
                  level =0.95, exact=FALSE)
test


