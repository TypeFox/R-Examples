### R code from vignette source 'laeken-variance.Rnw'

###################################################
### code chunk number 1: laeken-variance.Rnw:54-55
###################################################
options(prompt="R> ")


###################################################
### code chunk number 2: laeken-variance.Rnw:95-97 (eval = FALSE)
###################################################
## vignette("laeken-standard")
## vignette("laeken-pareto")


###################################################
### code chunk number 3: laeken-variance.Rnw:113-115
###################################################
library("laeken")
data("eusilc")


###################################################
### code chunk number 4: laeken-variance.Rnw:142-143
###################################################
args(variance)


###################################################
### code chunk number 5: laeken-variance.Rnw:228-231
###################################################
a <- arpr("eqIncome", weights = "rb050", data = eusilc)
variance("eqIncome", weights = "rb050", design = "db040", 
    data = eusilc, indicator = a, bootType = "naive", seed = 123)


###################################################
### code chunk number 6: laeken-variance.Rnw:239-242
###################################################
b <- arpr("eqIncome", weights = "rb050", breakdown = "db040", data = eusilc)
variance("eqIncome", weights = "rb050", breakdown = "db040", design = "db040", 
    data = eusilc, indicator = b, bootType = "naive", seed = 123)


###################################################
### code chunk number 7: laeken-variance.Rnw:311-314
###################################################
variance("eqIncome", weights = "rb050", design = "db040", 
    data = eusilc, indicator = a, X = calibVars(eusilc$db040), 
    seed = 123)


