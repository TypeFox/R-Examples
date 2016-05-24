### R code from vignette source 'equateIRT.Rnw'

###################################################
### code chunk number 1: equateIRT.Rnw:465-466
###################################################
options(prompt = "R> ", continue = "+ ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: equateIRT.Rnw:604-606
###################################################
library("equateIRT")
data("est2pl", package = "equateIRT")


###################################################
### code chunk number 3: equateIRT.Rnw:648-652
###################################################
test <- paste("test", 1:5, sep = "")
mod2pl <- modIRT(coef = est2pl$coef, var = est2pl$var, names = test, 
display = FALSE)
coef(mod2pl$test1)[1:5]


###################################################
### code chunk number 4: equateIRT.Rnw:662-663
###################################################
linkp(coef = est2pl$coef)


###################################################
### code chunk number 5: equateIRT.Rnw:679-681
###################################################
l15 <- direc(mod1 = mod2pl[1], mod2 = mod2pl[5], method = "Haebara")
l15


###################################################
### code chunk number 6: equateIRT.Rnw:692-693
###################################################
summary(l15)


###################################################
### code chunk number 7: equateIRT.Rnw:698-700
###################################################
direclist2pl <- alldirec(mods = mod2pl, method = "mean-mean")
direclist2pl


###################################################
### code chunk number 8: equateIRT.Rnw:704-705 (eval = FALSE)
###################################################
## summary(direclist2pl)


###################################################
### code chunk number 9: equateIRT.Rnw:708-709
###################################################
summary(direclist2pl, "test1.test5")


###################################################
### code chunk number 10: equateIRT.Rnw:718-720
###################################################
cec4 <- chainec(r = 4, direclist = direclist2pl)
cec4


###################################################
### code chunk number 11: equateIRT.Rnw:723-724 (eval = FALSE)
###################################################
## summary(cec4)


###################################################
### code chunk number 12: equateIRT.Rnw:728-729
###################################################
summary(cec4, "test1.test2.test3.test4")


###################################################
### code chunk number 13: equateIRT.Rnw:734-735
###################################################
cec4.1 <- chainec(r = 4, direclist = direclist2pl, f1 = "test1")


###################################################
### code chunk number 14: equateIRT.Rnw:741-743
###################################################
cec1234 <- chainec(r = 4, direclist = direclist2pl, f1 = "test1", 
f2 = "test4")


###################################################
### code chunk number 15: equateIRT.Rnw:748-752
###################################################
pth1 <- paste("test", c(1, 5, 4), sep = "")
pth1 <- data.frame(t(pth1), stringsAsFactors = FALSE)
cec154 <- chainec(direclist = direclist2pl, pths = pth1)
summary(cec154)


###################################################
### code chunk number 16: equateIRT.Rnw:756-760
###################################################
pth2 <- paste("test", 1:5, sep = "")
pth2 <- data.frame(t(pth2), stringsAsFactors = FALSE)
cec12345 <- chainec(direclist = direclist2pl, pths = pth2)
summary(cec12345)


###################################################
### code chunk number 17: equateIRT.Rnw:772-776
###################################################
ecall <- c(cec1234, cec154, cec12345, direclist2pl["test1.test5"])
fec <- bisectorec(ecall = ecall, mod = mod2pl, weighted = TRUE, 
unweighted = TRUE)
fec


###################################################
### code chunk number 18: equateIRT.Rnw:779-780
###################################################
summary(fec)


###################################################
### code chunk number 19: equateIRT.Rnw:786-787
###################################################
eqc(l15)


###################################################
### code chunk number 20: equateIRT.Rnw:790-791
###################################################
eqc(direclist2pl)


###################################################
### code chunk number 21: equateIRT.Rnw:794-795
###################################################
eqc(direclist2pl, link = "test1.test5")


###################################################
### code chunk number 22: equateIRT.Rnw:800-802
###################################################
eqc(cec4)
eqc(cec4, path = "test1.test2.test3.test4")


###################################################
### code chunk number 23: equateIRT.Rnw:807-809
###################################################
eqc(fec)
eqc(fec, link = "test1.test4", path = "bisector")


###################################################
### code chunk number 24: equateIRT.Rnw:817-820
###################################################
itm(l15)[1:3, ]
itm(direclist2pl, "test1.test5")[1:3, ]
itm(cec12345, "test1.test2.test3.test4.test5")[1:3, ]


###################################################
### code chunk number 25: equateIRT.Rnw:835-838
###################################################
eqc14 <- eqc(fec, link = "test1.test4", path = "bisector")
convert(A = eqc14$A, B = eqc14$B, coef = coef(mod2pl$test1), 
person.par = seq(-3, 3, 0.5))


###################################################
### code chunk number 26: equateIRT.Rnw:886-894
###################################################
library("ltm")
library("equateIRT")
data("data2pl", package = "equateIRT")
m1 <- ltm(data2pl[[1]] ~ z1)
m2 <- ltm(data2pl[[2]] ~ z1)
m3 <- ltm(data2pl[[3]] ~ z1)
m4 <- ltm(data2pl[[4]] ~ z1)
m5 <- ltm(data2pl[[5]] ~ z1)


###################################################
### code chunk number 27: equateIRT.Rnw:898-905
###################################################
estm1 <- import.ltm(m1, display = FALSE)
estm2 <- import.ltm(m2, display = FALSE)
estm3 <- import.ltm(m3, display = FALSE)
estm4 <- import.ltm(m4, display = FALSE)
estm5 <- import.ltm(m5, display = FALSE)
estm1$coef[1:3, ]
estm1$var[1:3,1:3]


###################################################
### code chunk number 28: equateIRT.Rnw:913-916
###################################################
estc <- list(estm1$coef, estm2$coef, estm3$coef, estm4$coef, estm5$coef)
estv <- list(estm1$var, estm2$var, estm3$var, estm4$var, estm5$var)
mod2pl.ltm <- modIRT(coef = estc, var = estv, display = FALSE)


