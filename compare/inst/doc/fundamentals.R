### R code from vignette source 'fundamentals.Rnw'

###################################################
### code chunk number 1: fundamentals.Rnw:17-19
###################################################
library(compare)



###################################################
### code chunk number 2: fundamentals.Rnw:37-39
###################################################
isTRUE(compare(1:10, 1:10))



###################################################
### code chunk number 3: fundamentals.Rnw:44-47
###################################################
obj1 <- c("a", "a", "b", "c")
obj1



###################################################
### code chunk number 4: fundamentals.Rnw:48-51
###################################################
obj2 <- factor(obj1)
obj2



###################################################
### code chunk number 5: fundamentals.Rnw:52-54
###################################################
transforms(compare(obj1, obj2[1:3], allowAll=TRUE))



###################################################
### code chunk number 6: fundamentals.Rnw:94-98
###################################################
compare(as.numeric(1:10), 
        as.character(10:1 + .1), 
        round=TRUE, coerce=TRUE, ignoreOrder=TRUE)



