### R code from vignette source 'chebpol.Rnw'

###################################################
### code chunk number 1: runge
###################################################
pdf('crunge.pdf')
library(chebpol)


###################################################
### code chunk number 2: runge
###################################################
f <- function(x) cos(3*pi*x)/(1+25*(x-0.25)^2)
ch <- Vectorize(chebappxf(f,15))
s <- seq(-1,1,length.out=401)
plot(s,f(s),type='l')
lines(s,ch(s), col='red')
kn <- chebknots(15)[[1]]
points(kn,f(kn))


###################################################
### code chunk number 3: runge
###################################################
invisible(dev.off())


###################################################
### code chunk number 4: runge
###################################################
pdf('ucrunge.pdf')
plot(s,f(s),type='l')
lines(s,ch(s), col='red')
points(kn,f(kn))


###################################################
### code chunk number 5: runge
###################################################
uc <- Vectorize(ucappxf(f,15))
lines(s,uc(s),col='blue')


###################################################
### code chunk number 6: runge
###################################################
invisible(dev.off())


###################################################
### code chunk number 7: multi
###################################################
library(chebpol)
f <- function(x) log(x[[1]])*sqrt(x[[2]])/log(sum(x))
ch <- chebappxf(f, c(5,8), list(c(1,2), c(15,20)))
uc <- ucappxf(f, c(5,8), list(c(1,2), c(15,20)))
tp <- c(runif(1,1,2), runif(1,15,20))
cat('arg:',tp,'true:', f(tp), 'ch:', ch(tp), 'uc:',uc(tp),'\n')


###################################################
### code chunk number 8: sqrt
###################################################
library(chebpol)
pdf('randgrid.pdf')


###################################################
### code chunk number 9: sqrt
###################################################
f <- function(x) cos(3*pi*x)/(1+25*(x-0.25)^2)
gr <- log(seq(exp(-1),exp(1),length=15))
chg <- Vectorize(chebappxgf(f,gr))
plot(s, f(s), col='black', type='l')
lines(s, chg(s), col='blue')
points(gr,f(gr))


###################################################
### code chunk number 10: runge
###################################################
invisible(dev.off())


###################################################
### code chunk number 11: ml
###################################################
pdf('ml.pdf')


###################################################
### code chunk number 12: ml
###################################################
f <- function(x) sign(sum(x^3)-0.1)*
                  sqrt(abs(25*prod(x)-4))/
                  (1+25*sum(x)^2)
grid <- replicate(4,list(seq(-1,1,length=15)))
ml <- mlappx(f,grid)
s <- seq(-1,1,length=400)
curve <- function(x) c(cos(1.2*pi*x),
                       sin(1.5*pi*x^3), 
                       x^2, -x/(1+x^2))
wf <- sapply(s,function(x) f(curve(x)))
wml <- sapply(s,function(x) ml(curve(x)))
plot(s,wf,typ='l')  # function
lines(s,wml,col='blue') # multilinear interpolation


###################################################
### code chunk number 13: ml
###################################################
invisible(dev.off())


###################################################
### code chunk number 14: poly
###################################################
set.seed(43)  # make sure we are reproducible


###################################################
### code chunk number 15: poly
###################################################
f <- function(x) 10/(10+sum(sqrt(x)))
knots <- matrix(runif(16000), 20)
phs <- polyh(f, knots, 3)
# pick some random point
x <- runif(20)
phs(x); f(x)


###################################################
### code chunk number 16: poly
###################################################
pdf('polyh.pdf')


###################################################
### code chunk number 17: poly
###################################################
f <- function(x) cos(3*pi*x)/(1+25*(x-0.25)^2)
knots <- runif(15, -1, 1)^3
phs <- polyh(f, knots, -36)
s <- seq(-1, 1, length=200)
plot(s, f(s), typ='l')
lines(s, phs(s), col='blue')
points(knots, f(knots))


###################################################
### code chunk number 18: poly
###################################################
invisible(dev.off())


