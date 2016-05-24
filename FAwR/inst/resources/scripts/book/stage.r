### R code from vignette source 'stage.rnw'

###################################################
### code chunk number 1: stage.rnw:5-7
###################################################
options(width=67)
source("../../scripts/functions.R")


###################################################
### code chunk number 2: get.stage
###################################################
stage <- read.csv("../../data/stage.csv")
str(stage)


###################################################
### code chunk number 3: stage.rnw:35-39
###################################################
stage$Tree.ID <- factor(stage$Tree.ID)
stage$Forest.ID <- factor(stage$Forest, labels = c("Kaniksu",
    "Coeur d'Alene", "St. Joe", "Clearwater", "Nez Perce", 
    "Clark Fork","Umatilla", "Wallowa", "Payette"))


###################################################
### code chunk number 4: stage.rnw:58-60
###################################################
stage$Hab.ID <- factor(stage$HabType, labels = c("Ts/Pac",
    "Ts/Op", "Th/Pach", "AG/Pach", "PA/Pach"))


###################################################
### code chunk number 5: stage.rnw:65-67
###################################################
stage$dbhib.cm <- stage$Dbhib * 2.54
stage$height.m <- stage$Height / 3.2808399


###################################################
### code chunk number 6: stage.rnw:73-74
###################################################
show.cols.with.na(stage)


