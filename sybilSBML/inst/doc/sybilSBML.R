### R code from vignette source 'sybilSBML.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: sybilSBML.Rnw:85-87 (eval = FALSE)
###################################################
## library(sybilSBML)
## model <- readSBMLmod("<model>.xml")


###################################################
### code chunk number 2: sybilSBML.Rnw:103-106
###################################################
library(sybilSBML)
mp     <- system.file(package = "sybilSBML", "extdata")
ec_mod <- file.path(mp, "ecoli_core_model.xml")


###################################################
### code chunk number 3: sybilSBML.Rnw:109-110
###################################################
mod <- readSBMLmod(ec_mod, bndCond = FALSE)


###################################################
### code chunk number 4: sybilSBML.Rnw:123-124
###################################################
err <- validateSBMLdocument(ec_mod)


