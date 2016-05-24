### R code from vignette source 'the_cwm_package.Rnw'

###################################################
### code chunk number 1: setup
###################################################
	options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
	set.seed(123)
	numSim=200


###################################################
### code chunk number 2: load
###################################################
library("cwm")


###################################################
### code chunk number 3: ex1
###################################################

library(MASS)
data(geyser)
x=geyser[,1]
y=geyser[,2]                            
cwrEmExample=cwrEm(x,y,nc=2)
print(cwrEmExample) 


