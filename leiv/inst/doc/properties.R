### R code from vignette source 'properties.Rnw'

###################################################
### code chunk number 1: data
###################################################
set.seed(1123)
n <- 20
X <- rnorm(n, mean=5, sd=4) # true x
x <- X + rnorm(n, mean=0, sd=2) # observed x
Y <- 2 + X # true y
y <- Y + rnorm(n, mean=0, sd=1) # observed y


###################################################
### code chunk number 2: original
###################################################
library(leiv)
l0 <- leiv(y ~ x)
print(l0)


###################################################
### code chunk number 3: suffstats
###################################################
r <- cor(x,y)
l <- sd(y)/sd(x)
xBar <- mean(x)
yBar <- mean(y)

l1 <- leiv(n=n, cor=r, sdRatio=l, xMean=xBar, yMean=yBar)
f <- function(beta) l1@density(beta)


###################################################
### code chunk number 4: ssplot
###################################################
plot(l0)
curve(f, add=TRUE, col="red", lty=2, lwd=2)


###################################################
### code chunk number 5: interchange
###################################################
l2 <- leiv(x ~ y) # slope is reciprocal slope of y vs. x
g <- function(beta) l2@density(1/beta)/beta^2


###################################################
### code chunk number 6: iplot
###################################################
plot(l0)
curve(g, add=TRUE, col="red", lty=2, lwd=2) 


###################################################
### code chunk number 7: scaling
###################################################
c <- 2
cy <- c*y
l3 <- leiv(cy ~ x) # slope is c times slope of y vs. x
h <- function(beta) l3@density(c*beta)*c


###################################################
### code chunk number 8: splot
###################################################
plot(l0)
curve(h, add=TRUE, col="red", lty=2, lwd=2) 


