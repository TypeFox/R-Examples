### R code from vignette source 'fia.rnw'

###################################################
### code chunk number 1: fia.rnw:4-5
###################################################
fia.plots <- read.csv("../../data/fia_plots.csv")


###################################################
### code chunk number 2: ratioSimulation
###################################################
fia.plots <- read.csv("../../data/fia_plots.csv")
fia.plots$forest <- factor(fia.plots$forest)
fia.plots$ba.m2.ha <- fia.plots$ba * 2.47105381 / 10.7639104
fia.plots$ht.m <- fia.plots$ht * 0.3048


