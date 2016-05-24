library(ic.infer)

mat <- as.matrix(swiss)
colnames(mat) <- NULL
covmat <- cov(swiss)
linmod <- lm(swiss)
linmodwt <- lm(swiss,weights=abs(-23:23))
linmodfac <- lm(1/time~poison+treat+poison:treat,boot::poisons)

## check normal functionality for covariance matrix
all.R2(covmat,ui=diag(1,5))
orlm1 <- orlm(covmat,ui=diag(c(-1,1,-1,1,1)),df.error=41)
summary(orlm1, overall.tests=FALSE)
ice <- ic.est(coef(linmod),Sigma=vcov(linmod),ui=(diag(c(-1,1,-1,1,1))),index=2:6)
summary(ice)
summary(ice,brief=TRUE)
orr1<-or.relimp(covmat,df.error=41,ui=diag(c(-1,1,-1,1,1)))

## check normal functionality for unnamed covariance matrix
covmat2 <- cov(mat)
all.R2(covmat2,ui=diag(1,5))
orlm1 <- orlm(covmat2,ui=diag(c(-1,1,-1,1,1)),df.error=41)
summary(orlm1, overall.tests=FALSE)
orr1<-or.relimp(covmat2,df.error=41,ui=diag(c(-1,1,-1,1,1)))

## check normal functionality for linear model
orlm1 <- orlm(linmod,ui=diag(c(-1,1,-1,1,1)))
summary(orlm1, overall.tests=FALSE)
orr1<-or.relimp(linmod,ui=diag(c(-1,1,-1,1,1)))
sum(orr1)
orlm1$R2

## check normal functionality for linear model of unnamed data
linmod2 <- lm(as.data.frame(mat))
orlm1 <- orlm(linmod2,ui=diag(c(-1,1,-1,1,1)))
summary(orlm1, overall.tests=FALSE)
orr1<-or.relimp(linmod2,df.error=41,ui=diag(c(-1,1,-1,1,1)))
sum(orr1)
orlm1$R2

## check normal functionality for weighted linear model
orlm1 <- orlm(linmodwt,ui=diag(c(-1,1,-1,1,1)))
summary(orlm1, overall.tests=FALSE)
orr1<-or.relimp(linmodwt,df.error=41,ui=diag(c(-1,1,-1,1,1)))
sum(orr1)
orlm1$R2

## check normal functionality for linear model with factors 
## and interactions
orlm1 <- orlm(linmodfac,ui=diag(c(1,1,1,-1,-1)),index=2:6)
orlm1
summary(orlm1, overall.tests=FALSE)
summary(orlm1,brief=TRUE, overall.tests=FALSE)
