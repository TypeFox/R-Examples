### R code from vignette source 'grade.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: intervals.ex1
###################################################
library(grade)
x <- c(1,2)
grade.interval(x, "[1,2]")
grade.interval(x, "[1,2.02]")
grade.interval(x, "[1,2.02]", tolerance=.03)

grade.interval("[NA, 1]", c(NA, 1), usena=T) # usena=T not allowed in grade.interval
grade.interval(c(1,Inf), c(1,Inf))
grade.interval(c(1,Inf), c(1,Inf), useinf=T)


###################################################
### code chunk number 2: <sets.ex1
###################################################
set1 <- "[1,2,3,5]"
set2 <- "[5,3,2,1]"
grade.orderedset(set1, set2)
grade.set(set1, set2)
set3 <- c(NA, Inf, pi)
grade.orderedset(set3, set3)
grade.orderedset(set3, set3, usena=T, useinf=T)


###################################################
### code chunk number 3: <discrete.ex1
###################################################
grade.discreteprobability(NULL, c(0, 1/2, 1/2, 0), checkcorrect=FALSE)
grade.discreteprobability(NULL, c(-1, 1/2, 1/2, 1), checkcorrect=FALSE)


###################################################
### code chunk number 4: <discrete.ex2
###################################################
correct <- c(0, 1/2, 1/4, 1/4)
sans <- "[0, 1/4, 1/4, 1/2]"
grade.discreteprobability(correct, sans)
grade.discreteprobability(correct, sans, ordered=T)


###################################################
### code chunk number 5: <neg.num.ex1
###################################################
grade.negative(NULL, -1)
grade.negative(NULL, -Inf)
grade.negative(NULL, -Inf, useinf=T)

grade.number(1, 1)
grade.number(3.141, "pi", tolerance=.002)


###################################################
### code chunk number 6: <parse.ex1
###################################################
grade.parse("[1, 2, 3]")
grade.parse("NA")
grade.parse("NA", usena=T)
grade.parse("[-Inf, Inf]")
grade.parse("[-Inf, Inf]", useinf=T)


###################################################
### code chunk number 7: <parse.ex2
###################################################
grade.parse(c(1,2,3))
grade.parse(NA)
grade.parse(c(-Inf, Inf), useinf=T)


###################################################
### code chunk number 8: <isscalar.ex1
###################################################
grade.isscalar(c(1,2))
grade.isscalar(Inf, useinf=T)


###################################################
### code chunk number 9: evalrm.ex1
###################################################
x <- NULL
grade.parse("x <- 1")
x
grade.parse("rm(list=ls())")
ls()
eval(parse(text="rm(list=ls())"))
ls()


###################################################
### code chunk number 10: <evalvsnum.ex1
###################################################
grade.parse("1")
grade.parse("1", useeval=F)
grade.parse("pi")
grade.parse("pi", useeval=F)
grade.parse("1/2")
grade.parse("1/2", useeval=F)


