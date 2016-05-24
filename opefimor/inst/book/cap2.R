###################################################
### chunk number 1: 
###################################################
#line 6 "cap2.Rnw"
options(prompt="R> ")
options(width=80)


###################################################
### chunk number 2: 
###################################################
#line 1120 "cap2.Rnw"
require(fBasics)

x <- seq(-10, 10, length=500)
y1 <- function(x)  dstable(x, alpha = 0.5, beta = 0.5, pm=1)
y2 <- function(x)  dstable(x, alpha = 0.75, beta = 0.5, pm=1)
y3 <- function(x)  dstable(x, alpha = 1, beta = 0.5, pm=1)
y4 <- function(x)  dstable(x, alpha = 1.25, beta = 0.5, pm=1)
y5 <- function(x)  dstable(x, alpha = 1.5, beta = 0.5, pm=1)

curve(y1, -5, 5, lty=6, ylim=c(0,0.6), ylab="density")
curve(y2, -5, 5, lty=2, add=TRUE)
curve(y3, -5, 5, lty=3, add=TRUE)
curve(y4, -5, 5, lty=4, add=TRUE)
curve(y5, -5, 5, lty=1, add=TRUE)

legend(-4, 0.6, legend=c(expression(alpha==0.5), 
expression(alpha==0.75), 
expression(alpha==1), 
expression(alpha==1.25), 
expression(alpha==1.5)), lty = c(6,2,3,4,1))


###################################################
### chunk number 3: 
###################################################
#line 1144 "cap2.Rnw"
require(fBasics)

x <- seq(-10, 10, length=500)
y1 <- function(x)  dstable(x, alpha = 0.5, beta = 0.5, pm=1)
y2 <- function(x)  dstable(x, alpha = 0.75, beta = 0.5, pm=1)
y3 <- function(x)  dstable(x, alpha = 1, beta = 0.5, pm=1)
y4 <- function(x)  dstable(x, alpha = 1.25, beta = 0.5, pm=1)
y5 <- function(x)  dstable(x, alpha = 1.5, beta = 0.5, pm=1)

curve(y1, -5, 5, lty=6, ylim=c(0,0.6), ylab="density")
curve(y2, -5, 5, lty=2, add=TRUE)
curve(y3, -5, 5, lty=3, add=TRUE)
curve(y4, -5, 5, lty=4, add=TRUE)
curve(y5, -5, 5, lty=1, add=TRUE)

legend(-4, 0.6, legend=c(expression(alpha==0.5), 
expression(alpha==0.75), 
expression(alpha==1), 
expression(alpha==1.25), 
expression(alpha==1.5)), lty = c(6,2,3,4,1))


###################################################
### chunk number 4: 
###################################################
#line 1173 "cap2.Rnw"
require(fBasics)

x <- seq(-10, 10, length=500)
y1 <- function(x)  dstable(x, alpha = 1, beta = -0.99, pm=1)
y2 <- function(x)  dstable(x, alpha = 1, beta = -0.3, pm=1)
y3 <- function(x)  dstable(x, alpha = 1, beta = 0, pm=1)
y4 <- function(x)  dstable(x, alpha = 1, beta = 0.5, pm=1)
y5 <- function(x)  dstable(x, alpha = 1, beta = 0.7, pm=1)

curve(y1, -3, 3, lty=6, ylim=c(0,0.35), ylab="density")
curve(y2, -3, 3, lty=2, add=TRUE)
curve(y3, -3, 3, lty=1, add=TRUE)
curve(y4, -3, 3, lty=4, add=TRUE)
curve(y5, -3, 3, lty=3, add=TRUE)

legend(-3, 0.33, legend=c(expression(beta==-0.99), 
expression(beta==-0.3), 
expression(beta==0), 
expression(beta==0.5), 
expression(beta==0.7)), lty = c(6,2,1,4,3))


###################################################
### chunk number 5: 
###################################################
#line 1198 "cap2.Rnw"
require(fBasics)

x <- seq(-10, 10, length=500)
y1 <- function(x)  dstable(x, alpha = 1, beta = -0.99, pm=1)
y2 <- function(x)  dstable(x, alpha = 1, beta = -0.3, pm=1)
y3 <- function(x)  dstable(x, alpha = 1, beta = 0, pm=1)
y4 <- function(x)  dstable(x, alpha = 1, beta = 0.5, pm=1)
y5 <- function(x)  dstable(x, alpha = 1, beta = 0.7, pm=1)

curve(y1, -3, 3, lty=6, ylim=c(0,0.35), ylab="density")
curve(y2, -3, 3, lty=2, add=TRUE)
curve(y3, -3, 3, lty=1, add=TRUE)
curve(y4, -3, 3, lty=4, add=TRUE)
curve(y5, -3, 3, lty=3, add=TRUE)


legend(-3, 0.33, legend=c(expression(beta==-0.99), 
expression(beta==-0.3), 
expression(beta==0), 
expression(beta==0.5), 
expression(beta==0.7)), lty = c(6,2,1,4,3))


###################################################
### chunk number 6: 
###################################################
#line 1294 "cap2.Rnw"
x <- seq(-3,3, length=20)
f <- function(x) dnorm(x)
f(x)


###################################################
### chunk number 7: 
###################################################
#line 1300 "cap2.Rnw"
y <- fft( f(x), inverse=TRUE)
y


###################################################
### chunk number 8: 
###################################################
#line 1306 "cap2.Rnw"
invFFT <- as.numeric( fft(y)/length(y) )
invFFT


###################################################
### chunk number 9: mle
###################################################
#line 2100 "cap2.Rnw"
set.seed(123)
library("stats4")
x <- rnorm(1000, mean=5, sd=2)
log.lik <- function(mu=1, sigma=1) 
 -sum( dnorm(x, mean=mu, sd=sigma,log=TRUE) )

fit <- mle(log.lik,lower=c(0,0),method="L-BFGS-B")
fit


###################################################
### chunk number 10: 
###################################################
#line 2111 "cap2.Rnw"
mean(x)
sd(x)


###################################################
### chunk number 11: 
###################################################
#line 2116 "cap2.Rnw"
logLik(fit)


###################################################
### chunk number 12: 
###################################################
#line 2120 "cap2.Rnw"
vcov(fit)


###################################################
### chunk number 13: 
###################################################
#line 2124 "cap2.Rnw"
confint(fit)
summary(fit)


