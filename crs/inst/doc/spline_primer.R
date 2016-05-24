### R code from vignette source 'spline_primer.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: spline_primer.Rnw:54-56
###################################################
library(crs)
options(prompt = "R> ", crs.messages = FALSE, digits = 3)


###################################################
### code chunk number 2: spline_primer.Rnw:133-139
###################################################
B <- function(x,b0,b1,b2) { (1-x)^2*b0+(1-x)*x*b1+x^2*b2 }
x <- seq(0,1,length=1000)
b0 <- 1
b1 <- -1
b2 <- 2
plot(x,B(x,b0,b1,b2),ylab="B(x)",type="l",lwd=2,cex.lab=1.25)


###################################################
### code chunk number 3: spline_primer.Rnw:216-221
###################################################
Bernstein <- function(n,i,x) { factorial(n)/(factorial(n-i)*factorial(i))*(1-x)^{n-i}*x^i }
x <- seq(0,1,length=100)
degree <- 2
plot(x,Bernstein(degree,0,x),type="l",lwd=2,ylab="B(x)",col=1)
for(i in 1:degree) lines(x,Bernstein(degree,i,x),lty=i+1,lwd=2,col=i+1)


###################################################
### code chunk number 4: spline_primer.Rnw:365-374
###################################################
degree <- 3
m <- degree+1
## nbreak is the total number of knots (2 indicates 2 endpoints, 0 interior)
nbreak <- 2+3
## N is the number of interior knots
N <- nbreak-2
x <- seq(0,1,length=1000)
B <- gsl.bs(x,degree=degree,nbreak=nbreak,intercept=TRUE)
matplot(x,B,type="l",lwd=2)


###################################################
### code chunk number 5: spline_primer.Rnw:376-379
###################################################
deriv <- 1
B.deriv <- gsl.bs(x,degree=degree,nbreak=nbreak,deriv=deriv,intercept=TRUE)
matplot(x,B.deriv,type="l",lwd=2)


