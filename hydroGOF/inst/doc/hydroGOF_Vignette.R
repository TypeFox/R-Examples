### R code from vignette source 'hydroGOF_Vignette.Rnw'

###################################################
### code chunk number 1: hydroGOF_Vignette.Rnw:22-23 (eval = FALSE)
###################################################
## install.packages("hydroGOF")


###################################################
### code chunk number 2: hydroGOF_Vignette.Rnw:33-34
###################################################
library(hydroGOF)


###################################################
### code chunk number 3: hydroGOF_Vignette.Rnw:39-42
###################################################
require(zoo)
data(EgaEnEstellaQts)
obs <- EgaEnEstellaQts


###################################################
### code chunk number 4: hydroGOF_Vignette.Rnw:47-48
###################################################
sim <- obs 


###################################################
### code chunk number 5: hydroGOF_Vignette.Rnw:53-54
###################################################
gof(sim=sim, obs=obs)


###################################################
### code chunk number 6: hydroGOF_Vignette.Rnw:59-60
###################################################
sim[1:2000] <- obs[1:2000] + rnorm(2000, mean=10)


###################################################
### code chunk number 7: hydroGOF_Vignette.Rnw:65-66
###################################################
ggof(sim=sim, obs=obs, ftype="dm", FUN=mean)


###################################################
### code chunk number 8: hydroGOF_Vignette.Rnw:78-79
###################################################
ggof(sim=sim, obs=obs, ftype="dm", FUN=mean, cal.ini="1963-01-01")


###################################################
### code chunk number 9: hydroGOF_Vignette.Rnw:84-88
###################################################
sim <- window(sim, start=as.Date("1963-01-01"))
obs <- window(obs, start=as.Date("1963-01-01"))

gof(sim, obs)


###################################################
### code chunk number 10: hydroGOF_Vignette.Rnw:104-105
###################################################
r <- sim-obs


###################################################
### code chunk number 11: hydroGOF_Vignette.Rnw:111-113
###################################################
library(hydroTSM)
smry(r) 


###################################################
### code chunk number 12: hydroGOF_Vignette.Rnw:116-118
###################################################
# daily, monthly and annual plots, boxplots and histograms
hydroplot(r, FUN=mean)


###################################################
### code chunk number 13: hydroGOF_Vignette.Rnw:123-125
###################################################
# daily, monthly and annual plots, boxplots and histograms
hydroplot(r, FUN=mean, pfreq="seasonal")


###################################################
### code chunk number 14: hydroGOF_Vignette.Rnw:136-139
###################################################
sessionInfo()$platform
sessionInfo()$R.version$version.string 
paste("hydroGOF", sessionInfo()$otherPkgs$hydroGOF$Version)


