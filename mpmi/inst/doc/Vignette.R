### R code from vignette source 'Vignette.Rnw'

###################################################
### code chunk number 1: Vignette.Rnw:18-19
###################################################
library(mpmi)


###################################################
### code chunk number 2: Vignette.Rnw:34-40
###################################################
# library(MASS)
# mu <- 1:100
# S <- toeplitz((100:1)/100)
# set.seed(123456789)
# dat <- mvrnorm(50, mu, S)
# cts <- scale(dat)


###################################################
### code chunk number 3: Vignette.Rnw:44-46
###################################################
data(mpmidata)
ctsresult <- cmi(cts)


###################################################
### code chunk number 4: Vignette.Rnw:51-52
###################################################
str(ctsresult)


###################################################
### code chunk number 5: Vignette.Rnw:55-56
###################################################
round(ctsresult$mi[1:5,1:5], 2)


###################################################
### code chunk number 6: Vignette.Rnw:59-60
###################################################
round(ctsresult$bcmi[1:5,1:5], 2)


###################################################
### code chunk number 7: Vignette.Rnw:66-67
###################################################
cmi.pw(cts[,1], cts[,1])


###################################################
### code chunk number 8: Vignette.Rnw:81-82
###################################################
mp(ctsresult$bcmi)


###################################################
### code chunk number 9: Vignette.Rnw:91-95
###################################################
# set.seed(987654321)
# disc <- rep(c("A", "H", "B"), ceiling(50 * 75 / 3))
# disc <- matrix(disc, nrow = 50, ncol = 75)
# disc <- apply(disc, 2, sample)


###################################################
### code chunk number 10: Vignette.Rnw:100-115
###################################################
cts2 <- cts
for (variable in 1:75)
{
    for (subject in 1:50)
    {
        if (disc[subject, variable] == "A") 
        {
            cts2[subject, variable] <- cts[subject, variable] - 2
        }
        if (disc[subject, variable] == "B") 
        {
            cts2[subject, variable] <- cts[subject, variable] - 2
        }
    }
}


###################################################
### code chunk number 11: Vignette.Rnw:118-119
###################################################
mixedresult <- mmi(cts2, disc)


###################################################
### code chunk number 12: Vignette.Rnw:127-128
###################################################
str(mixedresult, width = 60, strict.width = "cut")


###################################################
### code chunk number 13: Vignette.Rnw:131-132
###################################################
round(mixedresult$mi[1:5,1:5], 2)


###################################################
### code chunk number 14: Vignette.Rnw:135-136
###################################################
round(mixedresult$bcmi[1:5,1:5], 2)


###################################################
### code chunk number 15: Vignette.Rnw:140-141
###################################################
mmi.pw(cts2[,1], disc[,1])


###################################################
### code chunk number 16: Vignette.Rnw:147-148
###################################################
mp(mixedresult$bcmi)


###################################################
### code chunk number 17: Vignette.Rnw:179-181
###################################################
# library(parallel) # Commented for portability
library(compiler)


###################################################
### code chunk number 18: Vignette.Rnw:188-189
###################################################
hs <- apply(cts2, 2, dpik, level = 3L, kernel = "epanech")


###################################################
### code chunk number 19: Vignette.Rnw:196-206
###################################################
fi <- function(i)
{
    bcmis <- rep(NaN, 100)
    for (j in 1:100)
    {
        bcmis[j] <- mmi.pw(cts2[,j], disc[,i], h = hs[j])$bcmi
    }
    return(bcmis)
}
fi <- cmpfun(fi)


###################################################
### code chunk number 20: Vignette.Rnw:216-217
###################################################
# parmmi <- mcmapply(fi, 1:75)


###################################################
### code chunk number 21: Vignette.Rnw:221-222
###################################################
# sum(abs(mixedresult$bcmi - parmmi))


###################################################
### code chunk number 22: Vignette.Rnw:228-229
###################################################
hs2 <- apply(cts, 2, dpik, level = 3L, kernel = "epanech")


###################################################
### code chunk number 23: Vignette.Rnw:235-245
###################################################
fi <- function(i)
{
    bcmis <- rep(NaN, 100)
    for (j in i:100)
    {
        bcmis[j] <- cmi.pw(cts[,i], cts[,j], h = hs2[c(i,j)])$bcmi
    }
    return(bcmis)
}
fi <- cmpfun(fi)


###################################################
### code chunk number 24: Vignette.Rnw:252-253
###################################################
# parcmi <- mcmapply(fi, 1:100)


###################################################
### code chunk number 25: Vignette.Rnw:265-267
###################################################
lt <- function(x) x[lower.tri(x, diag = TRUE)]
# sum(abs(lt(ctsresult$bcmi) - lt(parcmi)))


