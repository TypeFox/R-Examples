### R code from vignette source 'kequate.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: kequate.Rnw:78-79
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: kequate.Rnw:414-417
###################################################
load("eqguide.RData")
load("CBsim.RData")
library(kequate)


###################################################
### code chunk number 3: kequate.Rnw:419-420
###################################################
freq <- kefreq(simeq$bivar1$X, 0:20)


###################################################
### code chunk number 4: kequate.Rnw:430-431
###################################################
SGfreq <- kefreq(simeq$bivar1$X, 0:20, simeq$bivar1$A, 0:10)


###################################################
### code chunk number 5: kequate.Rnw:433-436
###################################################
SGfreq <- kefreq(simeq$bivar1$X, 0:20, simeq$bivar1$A, 0:10)
PNEAT <- kefreq(simeq$bivar1$X, 0:20, simeq$bivar1$A, 0:10)
QNEAT <- kefreq(simeq$bivar2$Y, 0:20, simeq$bivar2$A, 0:10)


###################################################
### code chunk number 6: kequate.Rnw:447-449
###################################################
EGX <- glm(freq~I(X) + I(X^2) + I(X^3) + I(X^4) + I(X^5), family = 
"poisson", data = FXEG, x = TRUE)


###################################################
### code chunk number 7: kequate.Rnw:451-453
###################################################
EGY <- glm(freq~I(Y) + I(Y^2) + I(Y^3) + I(Y^4) + I(Y^5), family = 
"poisson", data = FYEG, x = TRUE)


###################################################
### code chunk number 8: kequate.Rnw:459-461
###################################################
SGglm <- glm(frequency~I(X) + I(X^2) + I(A) + I(A^2) + I(A^3) + I(X):I(A) 
+ I(X^2):I(A^2), data = SGfreq, family = "poisson", x = TRUE)


###################################################
### code chunk number 9: kequate.Rnw:465-469
###################################################
freqCB12 <- kefreq(CBeq12[,1], 0:40, CBeq12[,2])
freqCB21 <- kefreq(CBeq21[,1], 0:40, CBeq21[,2])
glmCB12 <- glm(frequency~I(X)+I(X^2)+I(X^3)+I(X^3)+I(Y)+I(Y^2)+I(Y^3)+I(Y^4)+I(X):I(Y)+I(X^2):I(Y)+I(X):I(Y^2)+I(X^2):I(Y^2), data=freqCB12, family=poisson, x=TRUE)
glmCB21 <- glm(frequency~I(X)+I(X^2)+I(X^3)+I(X^3)+I(Y)+I(Y^2)+I(Y^3)+I(Y^4)+I(X):I(Y)+I(X^2):I(Y)+I(X):I(Y^2)+I(X^2):I(Y^2), data=freqCB21, family=poisson, x=TRUE)


###################################################
### code chunk number 10: kequate.Rnw:489-490
###################################################
PNEATordered <- PNEAT[order(PNEAT$A, PNEAT$X),]


###################################################
### code chunk number 11: kequate.Rnw:493-496
###################################################
PNEAT$indx0 <- numeric(length(PNEAT$X))
PNEAT$ind1x <- numeric(length(PNEAT$X))
PNEAT$ind2x <- numeric(length(PNEAT$X))


###################################################
### code chunk number 12: kequate.Rnw:500-506
###################################################
PNEAT$indx0[PNEAT$X==0] <- 1
PNEAT$ind1x[PNEAT$X %in% c(5, 10, 15, 20)] <- 1
PNEAT$ind2x[PNEAT$X==5] <- 5
PNEAT$ind2x[PNEAT$X==10] <- 10
PNEAT$ind2x[PNEAT$X==15] <- 15
PNEAT$ind2x[PNEAT$X==20] <- 20


###################################################
### code chunk number 13: kequate.Rnw:508-518
###################################################
QNEAT$indy0 <- numeric(length(QNEAT$X))
QNEAT$ind1y <- numeric(length(QNEAT$X))
QNEAT$ind2y <- numeric(length(QNEAT$X))

QNEAT$indy0[QNEAT$X==0] <- 1
QNEAT$ind1y[QNEAT$X %in% c(5, 10, 15, 20)] <- 1
QNEAT$ind2y[QNEAT$X==5] <- 5
QNEAT$ind2y[QNEAT$X==10] <- 10
QNEAT$ind2y[QNEAT$X==15] <- 15
QNEAT$ind2y[QNEAT$X==20] <- 20


###################################################
### code chunk number 14: kequate.Rnw:522-525
###################################################
PNEATglm <- glm(frequency~I(X) + I(X^2) + I(X^3) + I(A) + I(A^2) + 
I(X):I(A) + I(X):I(A^2) + I(indx0) + I(ind1x) + I(ind2x) + I(ind2x^2), 
data = PNEAT, family = "poisson", x = TRUE)


###################################################
### code chunk number 15: kequate.Rnw:527-530
###################################################
QNEATglm <- glm(frequency~I(X) + I(X^2) + I(X^3) + I(A) + I(A^2) + 
I(X):I(A) + I(X):I(A^2) + I(indy0) + I(ind1y) + I(ind2y) + I(ind2y^2), 
data = QNEAT, family = "poisson", x = TRUE)


###################################################
### code chunk number 16: kequate.Rnw:536-537
###################################################
obs11 <- data11


###################################################
### code chunk number 17: kequate.Rnw:539-543 (eval = FALSE)
###################################################
## testfreq <- as.data.frame(table(factor(obs11$S11, levels = 0:40, ordered
## = TRUE), factor(obs11$edu, levels = 1:3, ordered = TRUE), 
## factor(obs11$math, levels = 1:3, ordered = TRUE), dnn = c("S11", "edu", 
## "math")))


###################################################
### code chunk number 18: kequate.Rnw:546-548 (eval = FALSE)
###################################################
## testdata11 <- data.frame(frequency = testfreq$Freq, S11 = rep(0:40, 6),
## edu = rep(1:2, each=41), math = rep(1:3, each = 41*2))


###################################################
### code chunk number 19: kequate.Rnw:550-552
###################################################
testdata11 <- data11
testdata12 <- data12


###################################################
### code chunk number 20: kequate.Rnw:556-559
###################################################
glm11 <- glm(frequency~I(S11) +  I(S11^2) + I(S11^3) + I(S11^4) + 
I(math) + I(math^2) + factor(edu) + I(S11):I(math) + I(S11):factor(edu) + 
I(math):factor(edu), data = testdata11, family = "poisson", x = TRUE)


###################################################
### code chunk number 21: kequate.Rnw:561-562
###################################################
glm12 <- glm(frequency~I(S12) +  I(S12^2) + I(S12^3) + I(S12^4) + I(math) + I(math^2) + factor(edu) + I(S12):I(math) + I(S12):factor(edu) + I(math):factor(edu), data = testdata12, family = "poisson", x = TRUE)


###################################################
### code chunk number 22: kequate.Rnw:571-572
###################################################
FTglm <- FTres(EGX$y, EGX$fitted.values)


###################################################
### code chunk number 23: kequate.Rnw:576-578
###################################################
Pest <- matrix(PNEATglm$fitted.values, nrow=21)
Pobs <- matrix(PNEATglm$y, nrow=21)


###################################################
### code chunk number 24: kequate.Rnw:580-581
###################################################
NEATPcdist <- cdist(Pest, Pobs)


###################################################
### code chunk number 25: kequate.Rnw:749-750
###################################################
keEG <- kequate("EG", 0:20, 0:20, EGX, EGY)


###################################################
### code chunk number 26: kequate.Rnw:754-756
###################################################
keEGobs <- kequate("EG", 0:20, 0:20, EGX$y/1453, EGY$y/1455, N = 1453, 
M = 1455, smoothed = FALSE)


###################################################
### code chunk number 27: kequate.Rnw:759-760
###################################################
summary(keEG)


###################################################
### code chunk number 28: kequate.Rnw:766-767
###################################################
keSG <- kequate("SG", 0:20, 0:10, SGglm)


###################################################
### code chunk number 29: kequate.Rnw:771-772
###################################################
summary(keSG)


###################################################
### code chunk number 30: kequate.Rnw:776-778
###################################################
DMSG <- SGglm$x[,-1]
PSG <- matrix(SGglm$fitted.values/sum(SGglm$fitted.values), nrow=21)


###################################################
### code chunk number 31: kequate.Rnw:780-781
###################################################
keSGDM <- kequate("SG", 0:20, 0:10, P = PSG, DM = DMSG, N = 1000)


###################################################
### code chunk number 32: kequate.Rnw:815-816
###################################################
keCB <- kequate("CB", 0:40, 0:40, glmCB12, glmCB21)


###################################################
### code chunk number 33: kequate.Rnw:821-829
###################################################
PNEAT <- kefreq(simeq$bivar1$X, 0:20, simeq$bivar1$A, 0:10)
QNEAT <- kefreq(simeq$bivar2$Y, 0:20, simeq$bivar2$A, 0:10)
NEATglmP <- glm(frequency~I(X) + I(X^2) + I(X^3) + I(A) + I(A^2) + 
I(X):I(A) + I(X):I(A^2), 
data = PNEAT, family = "poisson", x = TRUE)
NEATglmQ <- glm(frequency~I(X) + I(X^2) + I(X^3) + I(A) + I(A^2) + 
I(X):I(A) + I(X):I(A^2), 
data = QNEAT, family = "poisson", x = TRUE)


###################################################
### code chunk number 34: kequate.Rnw:833-834
###################################################
keNEATCE <- kequate("NEAT_CE", 0:20, 0:20, 0:10, NEATglmP, NEATglmQ)


###################################################
### code chunk number 35: neatceplot
###################################################
plot(keNEATCE)


###################################################
### code chunk number 36: neatceplot1
###################################################
plot(keNEATCE)


###################################################
### code chunk number 37: kequate.Rnw:853-854
###################################################
keNEATPSE <- kequate("NEAT_PSE", 0:20, 0:20, NEATglmP, NEATglmQ)


###################################################
### code chunk number 38: neatpseplot
###################################################
plot(keNEATPSE)


###################################################
### code chunk number 39: neatpseplot1
###################################################
plot(keNEATPSE)


###################################################
### code chunk number 40: kequate.Rnw:881-883
###################################################
keNEATPSEnew <- kequate("NEAT_PSE", 0:20, 0:20, PNEATglm, QNEATglm, hx = 
0.5, hy = 0.5, hxlin = 1000, hylin = 1000)


###################################################
### code chunk number 41: kequate.Rnw:887-888
###################################################
NECtest2012 <- kequate("NEC", 0:40, 0:40, glm12, glm11)


###################################################
### code chunk number 42: kequate.Rnw:891-892
###################################################
summary(NECtest2012)


###################################################
### code chunk number 43: necplot1
###################################################
plot(NECtest2012)


###################################################
### code chunk number 44: necplot
###################################################
plot(NECtest2012)


###################################################
### code chunk number 45: kequate.Rnw:908-910
###################################################
NECtestL <- kequate("NEC", 0:40, 0:40, glm12, glm11, kernel = "logistic")
NECtestU <- kequate("NEC", 0:40, 0:40, glm12, glm11, kernel = "uniform")


###################################################
### code chunk number 46: neccomp
###################################################
plot(0:40, getSee(NECtest2012), ylim=c(0, 0.8), pch=1, xlab="", ylab="")
par(new=TRUE)
plot(0:40, getSee(NECtestL), ylim=c(0, 0.8), pch=2, xlab="", ylab="")
par(new=TRUE)
plot(0:40, getSee(NECtestU), ylim=c(0, 0.8), pch=3, xlab="Score value", ylab="SEE")
legend("topright", inset=.1, title="Kernel utilized", c("Gaussian", "Logistic", "Uniform"), pch=c(1, 2, 3))


###################################################
### code chunk number 47: neckernelcomp
###################################################
plot(0:40, getSee(NECtest2012), ylim=c(0, 0.8), pch=1, xlab="", ylab="")
par(new=TRUE)
plot(0:40, getSee(NECtestL), ylim=c(0, 0.8), pch=2, xlab="", ylab="")
par(new=TRUE)
plot(0:40, getSee(NECtestU), ylim=c(0, 0.8), pch=3, xlab="Score value", ylab="SEE")
legend("topright", inset=.1, title="Kernel utilized", c("Gaussian", "Logistic", "Uniform"), pch=c(1, 2, 3))


###################################################
### code chunk number 48: kequate.Rnw:933-935
###################################################
keEGirt <- kequate("EG", 0:20, 0:20, EGX, EGY, irtx = simeq$irt2, irty = 
simeq$irt1)


###################################################
### code chunk number 49: seedplot
###################################################
SEEDPSECE <- genseed(keNEATPSE, keNEATCE)
plot(SEEDPSECE)


###################################################
### code chunk number 50: seedplot1
###################################################
SEEDPSECE <- genseed(keNEATPSE, keNEATCE)
plot(SEEDPSECE)


