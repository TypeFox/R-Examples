set.seed(1001)
library(regress)
n <- 101
x1 <- runif(n)
x2 <- seq(0,1,l=n)
z1 <- gl(4,10,n)
z2 <- gl(6,1,n)

X <- model.matrix(~1 + x1 + x2)
Z1 <- model.matrix(~z1-1)
Z2 <- model.matrix(~z2-1)

## Create the individual random and fixed effects
beta <- c(1,2,3)
eta1 <- rnorm(ncol(Z1),0,10)
eta2 <- rnorm(ncol(Z2),0,10)
eps <- rnorm(n,0,3)

## Combine them into a response
y <- X %*% beta + Z1 %*% eta1 + Z2 %*% eta2 + eps

## Fit the same model again
model <- regress(y~1 + x1 + x2,~z1 + z2, identity=TRUE,verbose=2)

summary(model)
names(model)

b1 <- BLUP(model,RE="z1")

b1$Mean
b1$Variance
b1$Covariance
cov2cor(b1$Covariance) ## Large correlation terms

## Simulate BLUPs given their mean and covariance structure
library(MASS)
tt <- mvrnorm(n=1000,mu=b1$Mean,Sigma=b1$Covariance)
plot(data.frame(tt))
matplot(t(tt),type="l",col="grey") ## So there is up and down shift only really.
lines(b1$Mean,col="red",lwd=2,type="b")

b2 <- BLUP(model,RE="z2")

b12 <- BLUP(model,RE=c("z1","z2")) ## To examine covariance of z1 and z2 BLUPs

bALL <- BLUP(model) ## To see blups of the measurement error epsilon

## Plot and compare the random effects and their BLUPs
plot(eta1,BLUP(model,RE="z1")$Mean,xlab="Random Effect",ylab="BLUP")
plot(eta2,BLUP(model,RE="z2")$Mean,xlab="Random Effect",ylab="BLUP")
plot(eps,BLUP(model,RE="In")$Mean,xlab="Random Effect",ylab="BLUP")

