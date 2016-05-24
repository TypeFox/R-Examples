### R code from vignette source 'qTableEx.Rnw'

###################################################
### code chunk number 1: qTableEx.Rnw:21-22
###################################################
options(continue = " ", digits = 5)


###################################################
### code chunk number 2: qTableEx.Rnw:34-40
###################################################
require("NMOF")
x <- rnorm(100L, mean = 0, sd = 1.5)
y <- rnorm(100L, mean = 1, sd = 1)
z <- rnorm(100L, mean = 1, sd = 0.5)
X <- cbind(x, y, z)
summary(X)


###################################################
### code chunk number 3: qTableEx.Rnw:44-46
###################################################
cat(qTable(X, yoffset = -0.025, unitlength = "5cm",
           circlesize = 0.0125, xmin = -10, xmax = 10, dec = 2))


###################################################
### code chunk number 4: res1
###################################################
## with limits
cat(qTable(X, yoffset = -0.025, unitlength = "5cm",
           circlesize = 0.0125, xmin = -10, xmax = 10, dec = 2))


###################################################
### code chunk number 5: res2
###################################################
## without specified limits
cat(qTable(X, yoffset = -0.025, unitlength = "5cm",
           circlesize = 0.0125, dec = 2))


###################################################
### code chunk number 6: res3
###################################################
## 3 digits
cat(qTable(X, yoffset = -0.025, unitlength = "5cm",
           circlesize = 0.0125, dec = 3))


###################################################
### code chunk number 7: res4
###################################################
## specific labels, but no limits
cat(qTable(X, yoffset = -0.025, unitlength = "5cm",
           labels = c(-8,2,8), at = c(-8,2,8),
           circlesize = 0.0125, dec = 1))


###################################################
### code chunk number 8: res5
###################################################
## specific labels and limits, linethickness
cat(qTable(X, yoffset = -0.025, unitlength = "5cm",
       labels = c("a","b","c"), at = c(-8,2,8),
       circlesize = 0.02, dec = 1, linethickness = "0.2ex",
       xmin = -10, xmax = 10))


###################################################
### code chunk number 9: res6
###################################################
## specific labels and limits, linethickness
cat(qTable(X, yoffset = -0.025, unitlength = "5cm",
       labels = c("a","b","c"), at = c(-8,2,8),
       circlesize = 0.02, dec = 1, linethickness = "0.2ex",
       xmin = -10, xmax = 10))


###################################################
### code chunk number 10: res7
###################################################
## with limits and alternative functions
cat(qTable(X, yoffset = -0.025, unitlength = "5cm",
           circlesize = 0.0125, xmin = -10, xmax = 10, dec = 2,
           funs = list(average = mean, 
                       `10th Q.` = function(x) quantile(x, 0.1),
                       `90th Q.` = function(x) quantile(x, 0.9))))


###################################################
### code chunk number 11: res8
###################################################
## with limits and without summary stats
cat(qTable(X, yoffset = -0.025, unitlength = "5cm",
           circlesize = 0.0125, xmin = -10, xmax = 10, dec = 2,
           funs = list()))


