library(ic.infer)

mat <- as.matrix(swiss)
colnames(mat) <- NULL
covmat <- cov(swiss)
linmod <- lm(swiss)
linmodwt <- lm(swiss,weights=abs(-23:23))
linmodfac <- lm(1/time~poison+treat+poison:treat,boot::poisons)

## test weights by simulation
## check normal functionality for covariance matrix
orlm1 <- orlm(covmat,ui=diag(c(-1,1,-1,1,1)),df.error=41)
summary(orlm1)
ic.test(orlm1)

## check normal functionality for unnamed covariance matrix
covmat2 <- cov(mat)
orlm1 <- orlm(covmat2,ui=diag(c(-1,1,-1,1,1)),df.error=41)
summary(orlm1)
ic.test(orlm1)

## additionally, bootstrap confidence intervals
## check normal functionality for linear model
orlm1 <- orlm(linmod,ui=diag(c(-1,1,-1,1,1)))
summary(orlm1)
ic.test(orlm1)
orlm1b <- orlm(linmod,ui=diag(c(-1,1,-1,1,1)),boot=TRUE,B=100)
summary(orlm1b)
ic.test(orlm1b)
orlm1bf <- orlm(linmod,ui=diag(c(-1,1,-1,1,1)),boot=TRUE,fixed=TRUE,B=100)
summary(orlm1bf)
ic.test(orlm1bf)

## check normal functionality for linear model of unnamed data
linmod2 <- lm(as.data.frame(mat))
orlm1 <- orlm(linmod2,ui=diag(c(-1,1,-1,1,1)))
summary(orlm1)
ic.test(orlm1)
orlm1b <- orlm(linmod2,ui=diag(c(-1,1,-1,1,1)),boot=TRUE,B=100)
summary(orlm1b)
ic.test(orlm1b)
orlm1bf <- orlm(linmod2,ui=diag(c(-1,1,-1,1,1)),boot=TRUE,fixed=TRUE,B=100)
summary(orlm1bf)
ic.test(orlm1bf)

## check normal functionality for weighted linear model
orlm1 <- orlm(linmodwt,ui=diag(c(-1,1,-1,1,1)))
summary(orlm1)
ic.test(orlm1)
orlm1b <- orlm(linmodwt,ui=diag(c(-1,1,-1,1,1)),boot=TRUE,B=100)
summary(orlm1b)
ic.test(orlm1b)
orlm1bf <- orlm(linmodwt,ui=diag(c(-1,1,-1,1,1)),boot=TRUE,fixed=TRUE,B=100)
summary(orlm1bf)
ic.test(orlm1bf)

## check normal functionality for linear model with factors 
## and interactions
orlm1 <- orlm(linmodfac,ui=diag(c(1,1,1,-1,-1)),index=2:6)
orlm1
summary(orlm1)
ic.test(orlm1)