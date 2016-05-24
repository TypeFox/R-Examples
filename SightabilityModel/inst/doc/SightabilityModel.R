### R code from vignette source 'SightabilityModel.Rnw'

###################################################
### code chunk number 1: SightabilityModel.Rnw:485-486
###################################################
 options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: SightabilityModel.Rnw:488-492
###################################################
 library("SightabilityModel")
 data("obs.m")
 data("exp.m")
 data("sampinfo.m")


###################################################
### code chunk number 3: SightabilityModel.Rnw:495-496
###################################################
exp.m[1:5, ] # first 5 observations


###################################################
### code chunk number 4: SightabilityModel.Rnw:506-507
###################################################
obs.m[1:5, ]


###################################################
### code chunk number 5: SightabilityModel.Rnw:524-525
###################################################
sampinfo.m


###################################################
### code chunk number 6: SightabilityModel.Rnw:535-538
###################################################
est.2004 <- Sight.Est(observed ~ voc, odat = subset(obs.m,
    year == 2004), sdat = exp.m, sampinfo = subset(sampinfo.m,
    year == 2004))


###################################################
### code chunk number 7: SightabilityModel.Rnw:547-548
###################################################
print(est.2004)


###################################################
### code chunk number 8: SightabilityModel.Rnw:567-568
###################################################
summary(est.2004)


###################################################
### code chunk number 9: SightabilityModel.Rnw:585-597
###################################################
tau.hats <- matrix(NA, 3, 5)
rownames(tau.hats) <- c("Stratum 1", "Stratum 2", "Stratum 3")
for(i in 1:3){
  tempsamp <- sampinfo.m[i, ]
  tempobs <- obs.m[obs.m$year == 2004 & obs.m$stratum == i, ]
  temp <- Sight.Est(observed ~ voc, odat = tempobs, sdat = exp.m,
      sampinfo = tempsamp)
  tau.hats[i, ] <- temp$est
}
colnames(tau.hats) <- names(temp$est)
tau.hats<-round(tau.hats, 0)
print(format(tau.hats, big.mark = ","),  quote = FALSE)


###################################################
### code chunk number 10: SightabilityModel.Rnw:602-609
###################################################
est.2004 <- Sight.Est(observed ~ voc, odat = subset(obs.m,
    year == 2004), sdat = exp.m, sampinfo = subset(sampinfo.m,
    year == 2004))
print(format(round(est.2004$est, 0), big.mark = ","), quote = FALSE)

naive.tau.hats<-round(apply(tau.hats[1:3, ], 2, sum), 0)
print(format(naive.tau.hats, big.mark = ","), quote = FALSE)


###################################################
### code chunk number 11: SightabilityModel.Rnw:642-651
###################################################
est.2006 <- Sight.Est(observed ~ voc, odat = subset(obs.m,
    year == 2006), sdat = exp.m, subset(sampinfo.m,
    year == 2006))
est.2007 <- Sight.Est(observed ~ voc, odat = subset(obs.m,
    year == 2007), sdat = exp.m, subset(sampinfo.m,
    year == 2007))

vdiff<-vardiff(est.2006, est.2007)
print(format(vdiff, nsmall = 0, big.mark = ","), quote = FALSE)


###################################################
### code chunk number 12: SightabilityModel.Rnw:654-656
###################################################
naive<-est.2006$est[2] + est.2007$est[2]
print(format(naive, nsmall = 0, big.mark = ","), quote = FALSE)


###################################################
### code chunk number 13: SightabilityModel.Rnw:688-689
###################################################
library("splines")


###################################################
### code chunk number 14: SightabilityModel.Rnw:692-699
###################################################
exp.m$voc.ns <- ns(exp.m$voc, df = 3)
obs.m$voc.ns <- predict(exp.m$voc.ns, obs.m$voc)
ns.est <- Sight.Est(observed ~ voc.ns, odat = subset(obs.m,
    year == 2004), sdat = exp.m, subset(sampinfo.m,
    year == 2004))
ns.est$sight
print(format(round(ns.est$est, 0),  big.mark = ","), quote = FALSE)


###################################################
### code chunk number 15: SightabilityModel.Rnw:726-727
###################################################
data("g.fit")


###################################################
### code chunk number 16: SightabilityModel.Rnw:729-730
###################################################
g.fit


###################################################
### code chunk number 17: SightabilityModel.Rnw:732-733
###################################################
data("gdat")


###################################################
### code chunk number 18: SightabilityModel.Rnw:735-736
###################################################
gdat[1:5, ]


###################################################
### code chunk number 19: SightabilityModel.Rnw:740-742
###################################################
sampinfo<-data.frame(nh = c(6, 23, 11), Nh =
c(6, 27, 65), stratum=c(1,2,3))


###################################################
### code chunk number 20: SightabilityModel.Rnw:746-747
###################################################
table(gdat$stratum)


###################################################
### code chunk number 21: SightabilityModel.Rnw:753-757
###################################################
goat.est<-Sight.Est(observed ~ GroupSize + Terrain + pct.VegCover,
    odat = gdat,  sampinfo = sampinfo[1:2, ], bet = g.fit$beta.g,
    varbet = g.fit$varbeta.g)
print(goat.est)


###################################################
### code chunk number 22: SightabilityModel.Rnw:785-796
###################################################
analytical.est <- Sight.Est(observed ~ voc, odat = subset(obs.m,
   year == 2004), sdat = subset(exp.m, year == 2005), subset(sampinfo.m,
   year == 2004), method = "Wong", logCI = T, alpha = 0.05,
   Vm.boot = FALSE)
print(format(round(analytical.est$est, 0), big.mark = ","), quote = FALSE)

boot.est <- Sight.Est(observed ~ voc, odat = subset(obs.m,
   year == 2004), sdat=subset(exp.m, year == 2005), subset(sampinfo.m,
   year == 2004), method = "Wong", logCI = T, alpha = 0.05,
   Vm.boot = TRUE, nboot = 10000)
print(format(round(boot.est$est, 0), big.mark = ","), quote = FALSE)


