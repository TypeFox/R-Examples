### R code from vignette source 'orsk_demo.Rnw'

###################################################
### code chunk number 1: orsk_demo.Rnw:157-158
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: orsk_demo.Rnw:160-161 (eval = FALSE)
###################################################
## install.packages("orsk")


###################################################
### code chunk number 3: orsk_demo.Rnw:164-166 (eval = FALSE)
###################################################
## library("orsk")
## vignette("orsk_demo",package = "orsk")


###################################################
### code chunk number 4: orsk_demo.Rnw:169-170 (eval = FALSE)
###################################################
## edit(vignette("orsk_demo",package = "orsk"))


###################################################
### code chunk number 5: orsk_demo.Rnw:180-181
###################################################
library("orsk")


###################################################
### code chunk number 6: orsk_demo.Rnw:183-186
###################################################
### analysis of Table 2 with grid method
res1 <- orsk(nctr=1636, ntrt=2601, a=2.61, al=2.25, au= 3.03, method="grid")
summary(res1)


###################################################
### code chunk number 7: orsk_demo.Rnw:193-194
###################################################
plot(res1, type="OR")


###################################################
### code chunk number 8: orsk_demo.Rnw:202-203
###################################################
plot(res1, type="RR")


###################################################
### code chunk number 9: orsk_demo.Rnw:222-226
###################################################
### analysis of Table 2 with optim method
require("setRNG")
old.seed <- setRNG(list(kind="Mersenne-Twister", normal.kind="Inversion",
seed=579))


###################################################
### code chunk number 10: orsk_demo.Rnw:228-229 (eval = FALSE)
###################################################
## res2 <- orsk(nctr=1636, ntrt=2601, a=2.61, al=2.25, au= 3.03, method="optim")


###################################################
### code chunk number 11: orsk_demo.Rnw:231-232
###################################################
res2 <- orsk(nctr=1636, ntrt=2601, a=2.61, al=2.25, au= 3.03, method="optim")


###################################################
### code chunk number 12: orsk_demo.Rnw:234-236
###################################################
summary(res2)
summary(res2$res$RR)


###################################################
### code chunk number 13: orsk_demo.Rnw:238-241
###################################################
#compare the computing speed between the two methods of estimation. Time to finish the modeling
#time.optim <- system.time(orsk(nctr=1636, ntrt=2601, a=2.61, al=2.25, au=3.03, method="optim"))[1]
time.grid <- system.time(orsk(nctr=1636, ntrt=2601, a=2.61, al=2.25, au=3.03, method="grid"))[1]


