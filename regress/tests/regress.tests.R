
## Test 1, Random Effects Model

library(nlme)
library(regress)
data(Oats)
names(Oats) <- c("B","V","N","Y")
Oats$N <- as.factor(Oats$N)

## Using regress
oats.reg <- regress(Y~N+V,~B+I(B:V),identity=TRUE,verbose=1,data=Oats)
summary(oats.reg)

## Using lme
oats.lme <- lme(Y~N+V,random=~1|B/V,data=Oats,method="REML")
summary(oats.lme)

if((oats.lme$sigma^2 - oats.reg$sigma[3])^2>0.0001) stop("Error 1 - doesn't match lme")

b1 <- ranef(oats.lme)
b2 <- BLUP(oats.reg)
if(sum((unlist(b1) - b2$Mean)^2)>0.0002) stop("Error with BLUP")

## Test 2, Multivariate Model

library(regress)
library(MASS)
n <- 100
mu <- c(1,2)
sigma <- c(10,10,5)
Sigma <- matrix(c(sigma[1],sigma[3],sigma[3],sigma[2]),2,2)
Y <- mvrnorm(n,mu,Sigma)

y <- c(Y[,1],Y[,2])
X <- kronecker(diag(1,2),rep(1,n))

V1 <- matrix(c(1,0,0,0),2,2)
V2 <- matrix(c(0,0,0,1),2,2)
V3 <- matrix(c(0,1,1,0),2,2)

sig1 <- kronecker(V1,diag(1,n))
sig2 <- kronecker(V2,diag(1,n))
gam <- kronecker(V3,diag(1,n))

reg.obj <- regress(y~X-1,~sig1+sig2+gam,identity=FALSE,verbose=1,start=c(10,10,5))
summary(reg.obj)
