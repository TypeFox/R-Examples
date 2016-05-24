### R code from vignette source 'correspondence.Rnw'

###################################################
### code chunk number 1: set_calculate_from_scratch
###################################################
calculate_from_scratch <- FALSE


###################################################
### code chunk number 2: useayl
###################################################
require(aylmer)
require(hyperdirichlet)


###################################################
### code chunk number 3: defa
###################################################
a <- matrix(c(26,1,2,5,18,9,4,4,NA),3,3,byrow=TRUE,dimnames=list("Reviewer 1"=c("yes","no", "missing"),"Reviewer 2"= c("yes", "no", "missing")))
aylmer_test_for_a <- aylmer.test(a,alternative=function(x){x[1,2]-x[2,1]})


###################################################
### code chunk number 4: aylmertesta
###################################################
a
aylmer.test(a,alternative=function(x)x[1,2]-x[2,1])


###################################################
### code chunk number 5: print.f1.aylmer.test
###################################################
f1 <- function(a)a[2,3]


###################################################
### code chunk number 6: correspondence.Rnw:194-195
###################################################
aylmer.test(a,alternative=f1)


###################################################
### code chunk number 7: calc.aaf1
###################################################
aaf1 <- aylmer.test(a,alternative=f1)


###################################################
### code chunk number 8: defb
###################################################
b <- uniform(4)
b[ 2] <-  19
b[ 3] <-   6
b[ 4] <-   9
b[ 5] <-   2
b[ 6] <-   4
b[ 9] <-  27
b[11] <-   4
b[13] <-   2


###################################################
### code chunk number 9: printb
###################################################
b


###################################################
### code chunk number 10: generate_data
###################################################
if(calculate_from_scratch){
  number <- 10000
  set.seed(0)
  d <- rhyperdirichlet(n=number, b)
  sam <- apply(d,1,function(x){x[2]/(x[2]+x[3])})
  rb <-  rbeta(number,2,6)
} else {
  load("precalc.Rdata")
}


###################################################
### code chunk number 11: sampat
###################################################
sampat <- apply(d,1,function(x){log(x[2]/x[3])})


###################################################
### code chunk number 12: plotPriors
###################################################

layout(matrix(1:4,2,2,byrow=TRUE))
hist(sam,freq=FALSE,col='gray',main='(a)',xlab=expression(psi))
abline(v=0.5,lwd=3,col='gray')
jj <- seq(from=0 , to=0.8 , len=100)
points(jj ,dbeta(jj,2,6),type='l')

par(pty='s')
qqplot(d,rb,asp=1,xlim=c(0,1),ylim=c(0,1),main='(b)',xlab=expression(paste(psi,' (all data)')),ylab=expression(paste(psi,' (complete cases)')))
abline(0,1)


par(pty='m')

jj.s <- seq_along(sam)/length(sam)
jj.r <- seq_along(rb)/length(rb)
  plot(sort(sam), jj.s, xlim=c(0,0.8),type='l',lty=1,main='(c)',xlab=expression(psi) , ylab='quantile')
points(sort(rb ), jj.r, type='l',lty=2)
legend('bottomright' , legend=c("c/cases", "all data"),lty=1:2)
segments(x0=0.5,y0=0.4,x1=0.5,y1=1, lwd=3,col='gray')

qqnorm(sampat,main='(d)',xlab='normal quantile' , ylab=expression(paste(phi, ' (empirical quantile)')))
abline(-1.819,var(sampat))


###################################################
### code chunk number 13: mlb
###################################################
mlb <- maximum_likelihood(b)


###################################################
### code chunk number 14: printmlb
###################################################
mlb


###################################################
### code chunk number 15: mlbf3
###################################################
f3 <- function(x){x[2]<x[3]}
mlb.f3 <- maximum_likelihood(b,disallowed=f3)


###################################################
### code chunk number 16: printmlb.f3
###################################################
mlb.f3


