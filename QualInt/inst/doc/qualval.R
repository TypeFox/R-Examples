### R code from vignette source 'qualval.Rnw'

###################################################
### code chunk number 1: qualval.Rnw:25-26
###################################################
library("QualInt")


###################################################
### code chunk number 2: qualval.Rnw:31-36
###################################################
test9 <- qualval(effect = c(1.0, 0.5, -2.0), 
                 se = c(0.86, 0.64, 0.32))
print(test9)
ibga(test9)
summary(test9)


###################################################
### code chunk number 3: qualval.Rnw:40-41
###################################################
plot(test9)


