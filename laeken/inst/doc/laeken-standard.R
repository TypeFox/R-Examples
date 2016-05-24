### R code from vignette source 'laeken-standard.Rnw'

###################################################
### code chunk number 1: laeken-standard.Rnw:54-55
###################################################
options(prompt="R> ")


###################################################
### code chunk number 2: laeken-standard.Rnw:137-139 (eval = FALSE)
###################################################
## vignette("laeken-pareto")
## vignette("laeken-variance")


###################################################
### code chunk number 3: laeken-standard.Rnw:152-155
###################################################
library("laeken")
data("eusilc")
head(eusilc, 3)


###################################################
### code chunk number 4: laeken-standard.Rnw:254-255
###################################################
methods(class="indicator")


###################################################
### code chunk number 5: laeken-standard.Rnw:333-335
###################################################
eusilc$eqSS <- eqSS("db030", "age", data=eusilc)
head(eusilc[,c("db030", "age", "eqSS")], 8)


###################################################
### code chunk number 6: laeken-standard.Rnw:347-354
###################################################
hplus <- c("hy040n", "hy050n", "hy070n", "hy080n", "hy090n", "hy110n")
hminus <- c("hy130n", "hy145n")
pplus <- c("py010n", "py050n", "py090n", "py100n", 
    "py110n", "py120n", "py130n", "py140n")
eusilc$eqIncome <- eqInc("db030", hplus, hminus, 
    pplus, character(), "eqSS", data=eusilc)
head(eusilc[,c("db030", "eqSS", "eqIncome")], 8)   


###################################################
### code chunk number 7: laeken-standard.Rnw:410-412
###################################################
weightedQuantile(eusilc$eqIncome, eusilc$rb050, 
    probs = c(0.2, 0.5, 0.8))


###################################################
### code chunk number 8: laeken-standard.Rnw:418-419
###################################################
weightedMedian(eusilc$eqIncome, eusilc$rb050)


###################################################
### code chunk number 9: laeken-standard.Rnw:431-433
###################################################
incMedian("eqIncome", weights = "rb050", data = eusilc) 
incQuintile("eqIncome", weights = "rb050", k = c(1, 4), data = eusilc)


###################################################
### code chunk number 10: laeken-standard.Rnw:523-525
###################################################
arpt("eqIncome", weights = "rb050", data = eusilc)
arpr("eqIncome", weights = "rb050", data = eusilc)


###################################################
### code chunk number 11: laeken-standard.Rnw:534-537
###################################################
arpr("eqIncome", weights = "rb050", p = 0.4, data = eusilc)
arpr("eqIncome", weights = "rb050", p = 0.5, data = eusilc)
arpr("eqIncome", weights = "rb050", p = 0.7, data = eusilc)


###################################################
### code chunk number 12: laeken-standard.Rnw:545-546
###################################################
arpr("eqIncome", weights = "rb050", breakdown = "db040", data = eusilc) 


###################################################
### code chunk number 13: laeken-standard.Rnw:554-557
###################################################
ageCat <- cut(eusilc$age, c(-1, 16, 25, 50, 65, Inf), right=FALSE)
eusilc$breakdown <- paste(ageCat, eusilc$rb090, sep=":")
arpr("eqIncome", weights = "rb050", breakdown = "breakdown", data = eusilc)     


###################################################
### code chunk number 14: laeken-standard.Rnw:593-594
###################################################
qsr("eqIncome", weights = "rb050", data = eusilc)   


###################################################
### code chunk number 15: laeken-standard.Rnw:600-601
###################################################
qsr("eqIncome", weights = "rb050", breakdown = "db040", data = eusilc) 


###################################################
### code chunk number 16: laeken-standard.Rnw:652-653
###################################################
rmpg("eqIncome", weights = "rb050", data = eusilc)


###################################################
### code chunk number 17: laeken-standard.Rnw:659-660
###################################################
rmpg("eqIncome", weights = "rb050", breakdown = "db040", data = eusilc)


###################################################
### code chunk number 18: laeken-standard.Rnw:667-670
###################################################
ageCat <- cut(eusilc$age, c(-1, 16, 25, 50, 65, Inf), right=FALSE)
eusilc$breakdown <- paste(ageCat, eusilc$rb090, sep=":")
rmpg("eqIncome", weights = "rb050", breakdown = "breakdown", data = eusilc)


###################################################
### code chunk number 19: laeken-standard.Rnw:694-695
###################################################
gini("eqIncome", weights = "rb050", data = eusilc)


###################################################
### code chunk number 20: laeken-standard.Rnw:700-701
###################################################
gini("eqIncome", weights = "rb050", breakdown = "db040", data = eusilc)


###################################################
### code chunk number 21: laeken-standard.Rnw:726-731
###################################################
a <- arpr("eqIncome", weights = "rb050", breakdown = "db040", data = eusilc)
print(a)
is.arpr(a)
is.indicator(a)
class(a)


###################################################
### code chunk number 22: laeken-standard.Rnw:741-742
###################################################
subset(a, strata = c("Lower Austria", "Vienna"))


