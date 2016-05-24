### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/PrcnApx.tex'

###################################################
### code chunk number 1: PrcnApx.tex:9-10
###################################################
library(HH)


###################################################
### code chunk number 2: PrcnApx.tex:13-14
###################################################
require(Rmpfr)


###################################################
### code chunk number 3: PrcnApx.tex:52-59
###################################################
## hhcapture("fp936.Rout", '
c(.9, (.3 + .6))
.9 == (.3 + .6)
.9 - (.3 + .6)
identical(.9, (.3 + .6))
all.equal(.9, (.3 + .6))
## ')


###################################################
### code chunk number 4: PrcnApx.tex:73-81
###################################################
## hhcapture("sqrt2.Rout", '
c(2, sqrt(2)^2)
sqrt(2)^2
2 == sqrt(2)^2
2 - sqrt(2)^2
identical(2, sqrt(2)^2)
all.equal(2, sqrt(2)^2)
## ')


###################################################
### code chunk number 5: PrcnApx.tex:182-185
###################################################
x <- 3.3125
print(x, digits=5)
sprintf("%+13.13a", x)


###################################################
### code chunk number 6: PrcnApx.tex:200-206
###################################################
## hhcapture("fpnums.Rout", '
nums <- c(0, .0625, .1, .3, .3125, .5, .6, (.3 + .6), .9,  1)
data.frame("decimal-2"=nums,
           "decimal-17"=format(nums, digits=17),
           hexadecimal=sprintf("%+13.13a", nums))
## ')


###################################################
### code chunk number 7: PrcnApx.tex:268-279
###################################################
old.width=options(width=80)
## hhcapture("as.bin.Rout", '
library(Rmpfr)
FourBits <- mpfr(matrix(0:39, 8, 5), precBits=4)
dimnames(FourBits) <- list(0:7, c(0,8,16,24,32))
FourBits
formatHex(FourBits)
formatBin(FourBits)
formatBin(FourBits, scientific=FALSE)
## ')
options(old.width)


###################################################
### code chunk number 8: PrcnApx.tex:374-383
###################################################
## hhcapture("fpnums369.Rout", '
nums369 <- c(.3, .6, .3+.6, 9)
nums369df <-
data.frame("decimal-2"=nums369,
           "decimal-17"=format(nums369, digits=17),
           hexadecimal=sprintf("%+13.13a", nums369))
nums369df[3,1] <- "0.3 + 0.6"
nums369df
## ')


###################################################
### code chunk number 9: PrcnApx.tex:414-417
###################################################
## hhcapture("sqrt2II.Rout", '
sprintf("%+13.13a", c(2, sqrt(2)^2))
## ')


###################################################
### code chunk number 10: PrcnApx.tex:438-442
###################################################
## hhcapture("machEps.Rout", '
c(100, 1e-10)
zapsmall(c(100, 1e-10))
## ')


###################################################
### code chunk number 11: PrcnApx.tex:469-507
###################################################
x <- 100
sprintf("%+13.13a", x)
sprintf("%+13.13a", (x+1))
sprintf("%+13.13a", (x+2))
sprintf("%+13.13a", x^2)
sprintf("%+13.13a", (x+1)^2)
sprintf("%+13.13a", (x+2)^2)
sprintf("%+13.13a", (x+2)^2 - (x+1)^2)
sprintf("%+13.13a", ((x+2)+(x+1)) * ((x+2)-(x+1)))

x
(x+1)
(x+2)
x^2
(x+1)^2
(x+2)^2
(x+2)^2 - (x+1)^2
((x+2)+(x+1)) * ((x+2)-(x+1))


x <- +0x1.p+27
sprintf("%+13.13a", x)
sprintf("%+13.13a", (x+1))
sprintf("%+13.13a", (x+2))
sprintf("%+13.13a", x^2)
sprintf("%+13.13a", (x+1)^2)
sprintf("%+13.13a", (x+2)^2)
sprintf("%+13.13a", (x+2)^2 - (x+1)^2)
sprintf("%+13.13a", ((x+2)+(x+1)) * ((x+2)-(x+1)))

x
(x+1)
(x+2)
x^2
(x+1)^2
(x+2)^2
(x+2)^2 - (x+1)^2
((x+2)+(x+1)) * ((x+2)-(x+1))


###################################################
### code chunk number 12: PrcnApx.tex:589-640
###################################################
## hhcapture("twopass.Rout", '
vartwo <- function(x) {
  n <- length(x)
  xbar <- mean(x)
  sum((x-xbar)^2) / (n-1)
}
x <- 1:3
vartwo(x)
vartwo(x+10^7)
## half machine precision
vartwo(x+10^8)
vartwo(x+10^15)
## boundary of machine precision
## See next table.
vartwo(x+10^16)
vartwo(x+10^17)
## ')

## hhcapture("onepass.Rout", '
varone <- function(x) {
  n <- length(x)
  xbar <- mean(x)
  (sum(x^2) - n*xbar^2) / (n-1)
}
x <- 1:3
varone(x)
varone(x+10^7)
## half machine precision
varone(x+10^8)
varone(x+10^15)
## boundary of machine precision
##
varone(x+10^16)
varone(x+10^17)
## ')

## hhcapture("twopassC.Rout", '
vartwoC <- function(x) {
  vartwo(x-mean(x))
}
x <- 1:3
vartwoC(x)
vartwoC(x+10^7)
## half machine precision
vartwoC(x+10^8)
vartwoC(x+10^15)
## boundary of machine precision
##
vartwoC(x+10^16)
vartwoC(x+10^17)
## ')


###################################################
### code chunk number 13: PrcnApx.tex:680-708
###################################################

## hhcapture("binvar53.Rout", '
x <- 1:3; p <- 15:16
xx <- t(outer(10^p, x, `+`)); dimnames(xx) <- list(x, p)
print(xx, digits=17)
formatHex(xx)

var(xx[,"15"])
var(xx[,"16"])

x54 <- mpfr(1:3, 54)
xx54 <- t(outer(10^p, x54, `+`)); dimnames(xx54) <- list(x, p)
xx54

vartwo(xx54[,"16"] - mean(xx54[,"16"]))

## hex for 54-bit numbers is not currently available from R.
## ')

hhcode("binvar53A.Rout", '
 > ## We manually constructed it here.
 > formatHex(xx54)  ## We manually constructed this
   15                     16
 1 +0x1.c6bf5263400080p+49 +0x1.1c37937e080008p+53
 2 +0x1.c6bf5263400100p+49 +0x1.1c37937e080010p+53
 3 +0x1.c6bf5263400180p+49 +0x1.1c37937e080018p+53
')



###################################################
### code chunk number 14: PrcnApx.tex:746-791
###################################################

## hhcapture("binvar5A.Rout", '
y <- 1:3; q <- 3:5; yy <- t(outer(2^q, y, `+`))
yy5 <- mpfr(yy, 5); dimnames(yy5) <- list(y, q)

yy5

vartwo(yy5[,"5"])

formatBin(yy5)

formatBin(yy5, scientific=FALSE)
## ')

## hhcapture("binvar5B.Rout", '
y <- 1:3; q <- 3:5; yy <- t(outer(2^q, y, `+`))
yy6 <- mpfr(yy, 6); dimnames(yy6) <- list(y, q)

yy6

vartwo(yy6[,"5"])

formatBin(yy6)

formatBin(yy6, scientific=FALSE)
## ')

## hhcapture("binvar5C.Rout", '
## the calculations in this chunk are not included in the book
yy5.num <- as.numeric(yy5)
yy5.num <- matrix(yy5.num, nrow(yy5), ncol(yy5), dimnames=dimnames(yy5))
apply(yy5.num, 2, var)

yy6.num <- as.numeric(yy6)
yy6.num <- matrix(yy6.num, nrow(yy6), ncol(yy6), dimnames=dimnames(yy6))
apply(yy6.num, 2, var)

vartwo(yy5[,1])
vartwo(yy5[,2])
vartwo(yy5[,3])

vartwo(yy6[,1])
vartwo(yy6[,2])
vartwo(yy6[,3])
## ')


###################################################
### code chunk number 15: PrcnApx.tex:878-897
###################################################
## hhcapture("binvar54.Rout", '
vartwo(xx54[,"15"])
## wrong answer.  numbers were shifted one binary position.
vartwo(xx54[,"16"])
## vartwoC protects against that problem and gets the right answer.
vartwoC(xx54[,"16"])

sum(xx54[1:2,"16"])

## Adding the first two numbers effectively doubled the numbers which
## means the significant bits were shifted one more place to the left.
## The first value was rounded up. Looking at just the last three bytes
## (where the last three bits are guaranteed 0):
##    +0x008p53 + 0x010p53 -> +0x010p53 + 0x010p53 -> +0x020p53

sum((xx54)[1:3,"16"])      ## too high
sum((xx54)[3:1,"16"])      ## too low
sum((xx54)[c(1,3,2),"16"]) ## just right
## ')


###################################################
### code chunk number 16: PrcnApx.tex:938-953
###################################################
## hhcapture("ModA.Rout", '
x <- 3; y <- 4

sqrt(x^2 + y^2)

Mod(x + 1i*y)

x <- 3e100; y <- 4e100
sqrt(x^2 + y^2)
Mod(x + y*1i)

x <- 3e305; y <- 4e305
sqrt(x^2 + y^2)
Mod(x + y*1i)
## ')


###################################################
### code chunk number 17: PrcnApx.tex:975-990
###################################################
## hhcapture("ModB.Rout", '
MyMod <- function(x, y) {
  XYmax <- max(abs(c(x, y)))
  xx <- x/XYmax
  yy <- y/XYmax

  result <- sqrt(xx^2 + yy^2)

  result * XYmax
}

x^2
y^2
MyMod(x, y)
## ')


###################################################
### code chunk number 18: PrcnApx.tex:1024-1050
###################################################
## hhcapture("onepassScalar.Rout", '
varoneScalar <- function(x) {
  ## This is a pedagogical example.
  ## Do not use this as a model for writing code.
  n <- length(x)
  sumx <- 0
  sumx2 <- 0
  for (i in 1:n) {
    sumx <- sumx + x[i]
    sumx2 <- sumx2 + x[i]^2
  }
    (sumx2 - (sumx^2)/n) / (n-1)
  }
x <- 1:3
varoneScalar(x)
varoneScalar(x+10^7)
## half machine precision
varoneScalar(x+10^8)

xx <- rnorm(1000)
## explicit loops are much slower in R
system.time(for (j in 1:1000) varoneScalar(xx))
system.time(for (j in 1:1000) varone(xx))
system.time(for (j in 1:1000) vartwo(xx))

## ')


