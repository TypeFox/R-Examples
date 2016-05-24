### R code from vignette source 'ShewhartConstants.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: ShewhartConstants.Rnw:11-30
###################################################
library(SixSigma)
nmax <- 25
n <- 2:nmax
d2 <- sapply(2:nmax, ss.cc.getd2)
d3 <- sapply(2:nmax, ss.cc.getd3)
c4 <- sapply(2:nmax, ss.cc.getc4)
A2 <- 3/(d2*sqrt(n))
D3 <- sapply(1:(nmax-1), function(x){
      max(c(0, 1 - 3*(d3[x]/d2[x])))})
D4 <- (1 + 3*(d3/d2))
B3 <- sapply(1:(nmax-1), function(x){
      max(0, 1 - 3*(sqrt(1-c4[x]^2)/c4[x]))})
B4 <- 1 + 3*(sqrt(1-c4^2)/c4)
constdf <- data.frame(n, d2, d3, c4, A2, 
    D3, D4, B3, B4)
print(xtable::xtable(constdf, 
    digits = 4,
    include.rownames = FALSE),
    floating = FALSE)    


