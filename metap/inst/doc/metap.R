### R code from vignette source 'metap.Rnw'

###################################################
### code chunk number 1: metap.Rnw:79-80
###################################################
library(metap)


###################################################
### code chunk number 2: metap.Rnw:94-98
###################################################
pvals <- c(0.1, 0.1, 0.9, 0.9, 0.9, 0.9)
istwo <- c(TRUE,  FALSE, TRUE, FALSE, TRUE, FALSE)
toinvert <- c(FALSE, TRUE, FALSE, FALSE, TRUE, TRUE)
two2one(pvals, two = istwo, invert = toinvert)


###################################################
### code chunk number 3: metap.Rnw:113-114
###################################################
data(validity)


###################################################
### code chunk number 4: metap.Rnw:118-119
###################################################
print(validity)


###################################################
### code chunk number 5: metap.Rnw:123-124
###################################################
schweder(validity)


###################################################
### code chunk number 6: metap.Rnw:148-151
###################################################
data(validity)
schweder(validity, drawline = c("bh", "ls", "ab"),
   ls.control = list(frac = 0.5), ab.control = list(a = 0, b = 0.01))


###################################################
### code chunk number 7: metap.Rnw:229-232
###################################################
sumlog(validity)
sumz(validity)
logitp(validity)


###################################################
### code chunk number 8: metap.Rnw:304-307
###################################################
minimump(validity)
sump(validity)
meanp(validity)


###################################################
### code chunk number 9: metap.Rnw:327-328
###################################################
votep(validity)


###################################################
### code chunk number 10: metap.Rnw:352-355
###################################################
pvals <- c(0.001, 0.001, 0.999, 0.999)
sumlog(pvals)
sumz(pvals)


