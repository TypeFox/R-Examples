#####################################
## Load stuff for testing purposes ##
#####################################

##- library(glmnet)
##- library(scalreg)
##- 
##- setwd("/u/meierluk/R/Pkgs/hdi/pkg/R")
##- 
##- source("hdi.R")
##- source("methods.R")
##- source("ridge-proj.R")
##- source("helpers.R")
##- 
##- ##########################
##- ## Load riboflavin data ##
##- ##########################
##- 
##- load("/u/meierluk/research/annualReview/Rcode/dsmN71.rda")

library(hdi)

data(riboflavin)

x <- riboflavin[,-1]
y <- riboflavin[,1]

dim(x)
##- [1]   71 4088
length(y)
##- [1] 71

x.use <- x[,1:200]

######################
## Ridge projection ##
######################

fit.ridge  <- ridge.proj(x = x.use, y = y, betainit = "scaled lasso")
fit.ridge1 <- ridge.proj(x = x.use + 10^-16, y = y, betainit = "scaled lasso")
stopifnot(all.equal(fit.ridge$pval, fit.ridge1$pval))

## Check standardization
fit.ridge2 <- ridge.proj(x = 2 + 4 * x.use, y = y, betainit = "scaled lasso")
stopifnot(all.equal(fit.ridge$pval, fit.ridge2$pval))

stopifnot(all.equal(max(abs(range(fit.ridge$bhat / fit.ridge2$bhat - 4))), 0))

## confidence intervals
ci.ridge <- confint(fit.ridge, level = 0.95)
ci.ridge2 <- confint(fit.ridge2, level = 0.95)

stopifnot(all.equal(ci.ridge, ci.ridge2 * 4))



#########
## GLM ##
#########



