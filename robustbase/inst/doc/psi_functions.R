### R code from vignette source 'psi_functions.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: init
###################################################
# set margins for plots
options(SweaveHooks=list(fig=function() par(mar=c(3,3,1.4,0.7),
                         mgp=c(1.5, 0.5, 0))))
require(robustbase)
source(system.file("xtraR/plot-psiFun.R", package = "robustbase", mustWork=TRUE))
##	       = ../xtraR/plot-psiFun.R
## x axis for plots:
x. <- seq(-5, 10, length=1501)


###################################################
### code chunk number 2: Huber
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(huberPsi, x., ylim=c(-1.4, 5), leg.loc="topright", main=FALSE)


###################################################
### code chunk number 3: lmrob-psi
###################################################
names(.Mpsi.tuning.defaults)


###################################################
### code chunk number 4: tuning-defaults
###################################################
print(c(k.M = .Mpsi.tuning.default("bisquare"),
        k.S = .Mchi.tuning.default("bisquare")), digits = 10)


###################################################
### code chunk number 5: bisquare
###################################################
getOption("SweaveHooks")[["fig"]]()
p.psiFun(x., "biweight", par = 4.685)


###################################################
### code chunk number 6: Hampel
###################################################
getOption("SweaveHooks")[["fig"]]()
## see also hampelPsi
p.psiFun(x., "Hampel", par = ## Default, but rounded:
                             round(c(1.5, 3.5, 8) * 0.9016085, 1))


###################################################
### code chunk number 7: GGW
###################################################
getOption("SweaveHooks")[["fig"]]()
p.psiFun(x., "GGW", par = c(-.5,1,.95,NA))


###################################################
### code chunk number 8: LQQ
###################################################
getOption("SweaveHooks")[["fig"]]()
p.psiFun(x., "LQQ", par = c(-.5,1.5,.95,NA))


###################################################
### code chunk number 9: optimal
###################################################
getOption("SweaveHooks")[["fig"]]()
p.psiFun(x., "optimal", par = 1.06, leg.loc="bottomright")


###################################################
### code chunk number 10: Welsh
###################################################
getOption("SweaveHooks")[["fig"]]()
p.psiFun(x., "Welsh", par = 2.11)


