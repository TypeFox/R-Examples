### R code from vignette source 'laeken-intro.Rnw'

###################################################
### code chunk number 1: laeken-intro.Rnw:108-110
###################################################
options(prompt = "R> ", continue = "+  ", width = 72, useFancyQuotes = FALSE)
library("laeken")


###################################################
### code chunk number 2: laeken-intro.Rnw:165-166 (eval = FALSE)
###################################################
## vignette(package="laeken")


###################################################
### code chunk number 3: laeken-intro.Rnw:245-247
###################################################
data("eusilc")
head(eusilc[, 1:10], 3)


###################################################
### code chunk number 4: laeken-intro.Rnw:270-272
###################################################
data("ses")
head(ses[, 1:7], 3)


###################################################
### code chunk number 5: laeken-intro.Rnw:393-394
###################################################
arpr("eqIncome", weights = "rb050", data = eusilc)


###################################################
### code chunk number 6: laeken-intro.Rnw:409-410
###################################################
arpr("eqIncome", weights = "rb050", p = c(0.4, 0.5, 0.7), data = eusilc)


###################################################
### code chunk number 7: laeken-intro.Rnw:432-433
###################################################
qsr("eqIncome", weights = "rb050", data = eusilc)


###################################################
### code chunk number 8: laeken-intro.Rnw:463-464
###################################################
rmpg("eqIncome", weights = "rb050", data = eusilc)


###################################################
### code chunk number 9: laeken-intro.Rnw:484-485
###################################################
gini("eqIncome", weights = "rb050", data = eusilc)


###################################################
### code chunk number 10: laeken-intro.Rnw:527-528
###################################################
gpg("earningsHour", gender = "sex", weigths = "weights", data = ses)


###################################################
### code chunk number 11: laeken-intro.Rnw:551-553
###################################################
gpg("earningsHour", gender = "sex", weigths = "weights", data = ses, 
    method = "median") 


###################################################
### code chunk number 12: laeken-intro.Rnw:590-591
###################################################
gini("eqIncome", weights = "rb050", data = eusilc)


###################################################
### code chunk number 13: laeken-intro.Rnw:594-595
###################################################
gini(eusilc$eqIncome, weights = eusilc$rb050)


###################################################
### code chunk number 14: laeken-intro.Rnw:671-673
###################################################
a <- arpr("eqIncome", weights = "rb050", breakdown = "db040", data = eusilc)
a


###################################################
### code chunk number 15: laeken-intro.Rnw:687-688
###################################################
subset(a, strata = c("Lower Austria", "Vienna"))


###################################################
### code chunk number 16: laeken-intro.Rnw:756-759
###################################################
hID <- eusilc$db030[which.max(eusilc$eqIncome)]
eqIncomeOut <- eusilc$eqIncome
eqIncomeOut[eusilc$db030 == hID] <- 10000000


###################################################
### code chunk number 17: laeken-intro.Rnw:766-768
###################################################
keep <- !duplicated(eusilc$db030)
eusilcH <- data.frame(eqIncome=eqIncomeOut, db090=eusilc$db090)[keep,]


###################################################
### code chunk number 18: laeken-intro.Rnw:797-798
###################################################
paretoQPlot(eusilcH$eqIncome, w = eusilcH$db090)


###################################################
### code chunk number 19: laeken-intro.Rnw:853-855
###################################################
ts <- paretoScale(eusilcH$eqIncome, w = eusilcH$db090)
ts


###################################################
### code chunk number 20: laeken-intro.Rnw:920-922
###################################################
thetaISE(eusilcH$eqIncome, k = ts$k, w = eusilcH$db090)
thetaISE(eusilcH$eqIncome, x0 = ts$x0, w = eusilcH$db090)


###################################################
### code chunk number 21: laeken-intro.Rnw:954-956
###################################################
thetaPDC(eusilcH$eqIncome, k = ts$k, w = eusilcH$db090)
thetaPDC(eusilcH$eqIncome, x0 = ts$x0, w = eusilcH$db090)


###################################################
### code chunk number 22: laeken-intro.Rnw:1010-1012
###################################################
fit <- paretoTail(eqIncomeOut, k = ts$k, w = eusilc$db090, 
                  groups = eusilc$db030)


###################################################
### code chunk number 23: laeken-intro.Rnw:1024-1025
###################################################
plot(fit)


###################################################
### code chunk number 24: laeken-intro.Rnw:1051-1053
###################################################
w <- reweightOut(fit, calibVars(eusilc$db040))
gini(eqIncomeOut, w)


###################################################
### code chunk number 25: laeken-intro.Rnw:1061-1064
###################################################
set.seed(123)
eqIncomeRN <- replaceOut(fit)
gini(eqIncomeRN, weights = eusilc$rb050)


###################################################
### code chunk number 26: laeken-intro.Rnw:1070-1072
###################################################
eqIncomeSN <- shrinkOut(fit)
gini(eqIncomeSN, weights = eusilc$rb050)


###################################################
### code chunk number 27: laeken-intro.Rnw:1079-1080
###################################################
gini(eqIncomeOut, weights = eusilc$rb050)


###################################################
### code chunk number 28: laeken-intro.Rnw:1153-1155
###################################################
arpr("eqIncome", weights = "rb050", design = "db040", cluster = "db030", 
     data = eusilc, var = "bootstrap", bootType = "naive", seed = 1234)


###################################################
### code chunk number 29: laeken-intro.Rnw:1203-1206
###################################################
aux <- cbind(calibVars(eusilc$db040), calibVars(eusilc$rb090))
arpr("eqIncome", weights = "rb050", design = "db040", cluster = "db030", 
     data = eusilc, var = "bootstrap", X = aux, seed = 1234)


