### R code from vignette source 'shortintro.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: check
###################################################
RNGkind()


###################################################
### code chunk number 2: example1
###################################################
RNGkind()
library(randtoolbox)
paramParkMiller <- c(mod=2^31-1, mult=16807, incr=0)
set.generator(name="congruRand", parameters=paramParkMiller, seed=1)
get.description()
RNGkind()
runif(10)


###################################################
### code chunk number 3: undo
###################################################
set.generator("default")
RNGkind()


###################################################
### code chunk number 4: example1
###################################################
setSeed(1)
congruRand(10, mod = 2^31-1, mult = 16807, incr = 0)


