### R code from vignette source 'agrmt.Rnw'

###################################################
### code chunk number 1: agrmt.Rnw:30-37
###################################################
library(agrmt)
x <- c(1,1,3) # these are our data
# 2 observations with position 1,
# 1 observation with position 3
table(x)
collapse(x)
collapse(x, pos=1:3) # now we specify which categories exist


###################################################
### code chunk number 2: agrmt.Rnw:85-95
###################################################
split.screen(c(2,2))
hist(c(1,1,1,1), breaks=0:4, xlim=c(0,4), ylim=c(0,4), col="grey", main="Agreement", xlab="Agreement = 1")
screen(2)
hist(c(1,2,3,4), breaks=0:4, xlim=c(0,4), ylim=c(0,4), col="grey", main="No Agreement", xlab="Agreement = 0")
screen(3)
hist(c(1,1,4,4), breaks=0:4, xlim=c(0,4), ylim=c(0,4), col="grey", main="Polarization", xlab="Agreement = -1")
screen(4)
hist(c(1,1,3,4), breaks=0:4, xlim=c(0,4), ylim=c(0,4), col="grey", main="", xlab="Agreement = 0.08")
# agreement(collapse(c(1,1,3,4), pos=1:4))
close.screen(all=TRUE)


###################################################
### code chunk number 3: agrmt.Rnw:104-106
###################################################
x <- c(30, 40, 210, 130, 530, 50, 10) # these are our data
agreement(x)


###################################################
### code chunk number 4: agrmt.Rnw:121-123
###################################################
agreement(c(2.4,2.8,3.2,6.2,13.5,30.4,41.6)) # PvdA
agreement(c(1.6,2.6,8.2,21,29.3,27,10.3))    # D66


###################################################
### code chunk number 5: agrmt.Rnw:146-148
###################################################
polarization(collapse(c(1,2,4,2,5,2,7,7,3,1,2,1,3,2,4,
1,5,2,3,2,4,2,3,1,1,3), pos=1:7))


###################################################
### code chunk number 6: agrmt.Rnw:153-155
###################################################
polarization(c(2.4,2.8,3.2,6.2,13.5,30.4,41.6)) # PvdA
polarization(c(1.6,2.6,8.2,21,29.3,27,10.3))    # D66


###################################################
### code chunk number 7: agrmt.Rnw:217-231
###################################################
# Example 1: different types of distributions
V1 <- c(0,1,0)     # A
ajus(V1)
V2 <- c(0,0,1)     # J
ajus(V2)
V3 <- c(1,0,1)     # U
ajus(V3)
V4 <- c(1,0,1,0)   # S
ajus(V4)
V5 <- c(0,0,0)     # F
ajus(V5)
V6 <- c(1,0,0)     # L
ajus(V6)
ajus(V6, variant="strict") # gives J


###################################################
### code chunk number 8: agrmt.Rnw:236-241
###################################################
# Example 2: varying tolerance to check sensitivity
V7 <- c(0,0,1,2,1)
ajus(V7, tolerance=0.5)
ajus(V7, tolerance=1)
ajusCheck(V7, t=c(0.1, 0.5, 0.6, 1))


###################################################
### code chunk number 9: agrmt.Rnw:244-246
###################################################
# Example 3: plotting AJUS
ajusPlot(V7)


###################################################
### code chunk number 10: agrmt.Rnw:278-281
###################################################
# Example 1: finding the mode
V1 <- c(30,40,210,130,530,50,10)
modes(V1) # will find position 5


###################################################
### code chunk number 11: agrmt.Rnw:286-289
###################################################
# Example 2:
V2 <- c(3,0,4,1)
modes(V2) # will find position 3


###################################################
### code chunk number 12: agrmt.Rnw:294-296
###################################################
# Example 3: providing categories
modes(V2,pos=-1:2) # will still find position 3, but give the value of 1 as mode


###################################################
### code chunk number 13: agrmt.Rnw:301-304
###################################################
# Example 4: similar values
V3 <- c(30,40,510,130,530,50,10)
modes(V3, tolerance=30) # will find positions 3 and 5 (510 and 530 are nearly the same)


