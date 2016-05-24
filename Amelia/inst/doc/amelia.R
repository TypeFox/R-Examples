### R code from vignette source 'amelia.Rnw'

###################################################
### code chunk number 1: amelia.Rnw:78-79
###################################################
ameliaVersion <- packageVersion("Amelia")


###################################################
### code chunk number 2: amelia.Rnw:492-496
###################################################
options("digits"=4)
options("width"=70)
options("show.signif.stars" = FALSE)
set.seed(12345)


###################################################
### code chunk number 3: amelia.Rnw:501-503
###################################################
require(Amelia)
data(freetrade)


###################################################
### code chunk number 4: amelia.Rnw:509-510
###################################################
summary(freetrade)


###################################################
### code chunk number 5: amelia.Rnw:517-519
###################################################
summary(lm(tariff ~ polity + pop + gdp.pc + year + country,
          data = freetrade))


###################################################
### code chunk number 6: amelia.Rnw:555-557
###################################################
a.out <- amelia(freetrade, m = 5, ts = "year", cs = "country")
a.out


###################################################
### code chunk number 7: hist1plot
###################################################
hist(a.out$imputations[[3]]$tariff, col="grey", border="white")


###################################################
### code chunk number 8: amelia.Rnw:579-580
###################################################
options(SweaveHooks = list(fig = function() par(mfrow=c(1,1))))


###################################################
### code chunk number 9: hist1
###################################################
getOption("SweaveHooks")[["fig"]]()
hist(a.out$imputations[[3]]$tariff, col="grey", border="white")


###################################################
### code chunk number 10: amelia.Rnw:652-654
###################################################
a.out.more <- amelia(freetrade, m = 10, ts = "year", cs = "country", p2s=0)
a.out.more


###################################################
### code chunk number 11: amelia.Rnw:658-660
###################################################
a.out.more <- ameliabind(a.out, a.out.more)
a.out.more


###################################################
### code chunk number 12: amelia.Rnw:712-713
###################################################
amelia(freetrade, m = 1, ts = "year", cs = "country", p2s = 2)


###################################################
### code chunk number 13: amelia.Rnw:813-814
###################################################
table(a.out$imputations[[3]]$polity)


###################################################
### code chunk number 14: amelia.Rnw:830-833
###################################################
a.out1 <- amelia(freetrade, m = 5, ts = "year", cs = "country", ords =
                 "polity", p2s = 0)
table(a.out1$imputations[[3]]$polity)


###################################################
### code chunk number 15: amelia.Rnw:850-851
###################################################
table(a.out1$imputations[[3]]$signed)


###################################################
### code chunk number 16: amelia.Rnw:868-871
###################################################
a.out2 <- amelia(freetrade, m = 5, ts = "year", cs = "country", noms =
                 "signed", p2s = 0)
table(a.out2$imputations[[3]]$signed)


###################################################
### code chunk number 17: logshist
###################################################
hist(freetrade$tariff, col="grey", border="white")
hist(log(freetrade$tariff), col="grey", border="white")


###################################################
### code chunk number 18: amelia.Rnw:905-906
###################################################
options(SweaveHooks = list(fig = function() par(mfrow=c(1,2))))


###################################################
### code chunk number 19: hist2
###################################################
getOption("SweaveHooks")[["fig"]]()
hist(freetrade$tariff, col="grey", border="white")
hist(log(freetrade$tariff), col="grey", border="white")


###################################################
### code chunk number 20: amelia.Rnw:920-921
###################################################
options(SweaveHooks = list(fig = function() par(mfrow=c(1,1))))


###################################################
### code chunk number 21: amelia.Rnw:956-957
###################################################
amelia(freetrade, idvars = c("year", "country"))


###################################################
### code chunk number 22: amelia.Rnw:964-965
###################################################
a.out2 <- amelia(freetrade, idvars = c("year"))


###################################################
### code chunk number 23: amelia.Rnw:1001-1002
###################################################
a.out2 <- amelia(freetrade, ts = "year", cs = "country", polytime = 2)


###################################################
### code chunk number 24: amelia.Rnw:1022-1024
###################################################
a.out.time <- amelia(freetrade, ts = "year", cs = "country", polytime = 2,
                 intercs = TRUE, p2s = 2)


###################################################
### code chunk number 25: tcomp1
###################################################
tscsPlot(a.out, cs = "Malaysia", main = "Malaysia (no time settings)",
         var = "tariff", ylim = c(-10, 60))

tscsPlot(a.out.time, cs = "Malaysia", main = "Malaysia (with time settings)",
         var = "tariff", ylim = c(-10, 60))


###################################################
### code chunk number 26: amelia.Rnw:1044-1045
###################################################
options(SweaveHooks = list(fig = function() par(mfrow=c(1,2))))


###################################################
### code chunk number 27: timecompare
###################################################
getOption("SweaveHooks")[["fig"]]()
tscsPlot(a.out, cs = "Malaysia", main = "Malaysia (no time settings)",
         var = "tariff", ylim = c(-10, 60))

tscsPlot(a.out.time, cs = "Malaysia", main = "Malaysia (with time settings)",
         var = "tariff", ylim = c(-10, 60))


###################################################
### code chunk number 28: amelia.Rnw:1077-1079
###################################################
a.out2 <- amelia(freetrade, ts = "year", cs = "country", lags = "tariff",
                 leads = "tariff")


###################################################
### code chunk number 29: amelia.Rnw:1106-1107
###################################################
a.out.time


###################################################
### code chunk number 30: amelia.Rnw:1130-1133
###################################################
a.out.time2 <- amelia(freetrade, ts = "year", cs = "country", polytime = 2,
                 intercs = TRUE, p2s = 0, empri = .01*nrow(freetrade))
a.out.time2


###################################################
### code chunk number 31: tcomp2
###################################################
tscsPlot(a.out.time, cs = "Malaysia", main = "Malaysia (no ridge prior)",
         var = "tariff", ylim = c(-10, 60))

tscsPlot(a.out.time2, cs = "Malaysia", main = "Malaysia (with ridge prior)",
         var = "tariff", ylim = c(-10, 60))


###################################################
### code chunk number 32: amelia.Rnw:1146-1147
###################################################
options(SweaveHooks = list(fig = function() par(mfrow=c(1,2))))


###################################################
### code chunk number 33: timecomp2
###################################################
getOption("SweaveHooks")[["fig"]]()
tscsPlot(a.out.time, cs = "Malaysia", main = "Malaysia (no ridge prior)",
         var = "tariff", ylim = c(-10, 60))

tscsPlot(a.out.time2, cs = "Malaysia", main = "Malaysia (with ridge prior)",
         var = "tariff", ylim = c(-10, 60))


###################################################
### code chunk number 34: amelia.Rnw:1201-1202
###################################################
freetrade[freetrade$country == "Thailand", c("year","country","tariff")]


###################################################
### code chunk number 35: amelia.Rnw:1211-1213
###################################################
pr <- matrix(c(158,159,160,3,3,3,40,40,40,3,3,3), nrow=3, ncol=4)
pr


###################################################
### code chunk number 36: amelia.Rnw:1220-1221
###################################################
a.out.pr <- amelia(freetrade, ts = "year", cs = "country", priors = pr)


###################################################
### code chunk number 37: amelia.Rnw:1230-1232
###################################################
pr.2 <- matrix(c(158,159,160,3,3,3,34,34,34,46,46,46,.95,.95,.95), nrow=3, ncol=5)
pr.2


###################################################
### code chunk number 38: amelia.Rnw:1243-1245
###################################################
pr.3 <- matrix(c(158,159,160,0,3,3,3,3,40,40,40,20,3,3,3,5), nrow=4, ncol=4)
pr.3


###################################################
### code chunk number 39: amelia.Rnw:1289-1291
###################################################
bds <- matrix(c(3, 30, 40), nrow = 1, ncol = 3)
bds


###################################################
### code chunk number 40: amelia.Rnw:1296-1298
###################################################
a.out.bds <- amelia(freetrade, ts = "year", cs = "country", bounds = bds,
                    max.resample = 1000)


###################################################
### code chunk number 41: bounds
###################################################
tscsPlot(a.out, cs = "Malaysia", main = "No logical bounds", var =
         "tariff", ylim = c(-10,60))

tscsPlot(a.out.bds, cs = "Malaysia", main = "Bounded between 30 and 40", var =
         "tariff", ylim = c(-10,60))


###################################################
### code chunk number 42: amelia.Rnw:1321-1322
###################################################
options(SweaveHooks = list(fig = function() par(mfrow=c(1,2))))


###################################################
### code chunk number 43: boundscomp
###################################################
getOption("SweaveHooks")[["fig"]]()
tscsPlot(a.out, cs = "Malaysia", main = "No logical bounds", var =
         "tariff", ylim = c(-10,60))

tscsPlot(a.out.bds, cs = "Malaysia", main = "Bounded between 30 and 40", var =
         "tariff", ylim = c(-10,60))


###################################################
### code chunk number 44: plotmeth
###################################################
plot(a.out, which.vars = 3:6)


###################################################
### code chunk number 45: plot1
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(a.out, which.vars = 3:6)


###################################################
### code chunk number 46: amelia.Rnw:1396-1397
###################################################
compare.density(a.out, var = "signed")


###################################################
### code chunk number 47: amelia.Rnw:1400-1401
###################################################
options(SweaveHooks = list(fig = function() par(mfrow=c(1,1))))


###################################################
### code chunk number 48: overimp
###################################################
overimpute(a.out, var = "tariff")


###################################################
### code chunk number 49: amelia.Rnw:1434-1435
###################################################
options(SweaveHooks = list(fig = function() par(mfrow=c(1,1))))


###################################################
### code chunk number 50: oi2
###################################################
getOption("SweaveHooks")[["fig"]]()
overimpute(a.out, var = "tariff")


###################################################
### code chunk number 51: overimp-bad
###################################################
dd <- Amelia:::rmvnorm(50, mu = c(0.5,0.5), vcv =
                       matrix(c(0.25^2,.06, .06,0.25^2),2,2))
ddmiss <- sample(1:50, replace = FALSE, size = 10)
is.na(dd) <- ddmiss
aa.out <- amelia(dd, m= 5)
overimpute(aa.out, var = 2, main = "Observed versus Imputed Values")


###################################################
### code chunk number 52: oi
###################################################
getOption("SweaveHooks")[["fig"]]()
dd <- Amelia:::rmvnorm(50, mu = c(0.5,0.5), vcv =
                       matrix(c(0.25^2,.06, .06,0.25^2),2,2))
ddmiss <- sample(1:50, replace = FALSE, size = 10)
is.na(dd) <- ddmiss
aa.out <- amelia(dd, m= 5)
overimpute(aa.out, var = 2, main = "Observed versus Imputed Values")


###################################################
### code chunk number 53: disp1d
###################################################
disperse(a.out, dims = 1, m = 5)
disperse(a.out, dims = 2, m = 5)


###################################################
### code chunk number 54: amelia.Rnw:1560-1561
###################################################
options(SweaveHooks = list(fig = function() par(mfrow=c(1,2))))


###################################################
### code chunk number 55: disp1dfig
###################################################
getOption("SweaveHooks")[["fig"]]()
disperse(a.out, dims = 1, m = 5)
disperse(a.out, dims = 2, m = 5)


###################################################
### code chunk number 56: amelia.Rnw:1649-1650
###################################################
options(SweaveHooks = list(fig = function() par(mfrow=c(1,1))))


###################################################
### code chunk number 57: tsplot1
###################################################
tscsPlot(a.out.time, cs = "Malaysia", main = "Malaysia (with time settings)",
         var = "tariff", ylim = c(-10, 60))


###################################################
### code chunk number 58: tsplot2
###################################################
getOption("SweaveHooks")[["fig"]]()
tscsPlot(a.out.time, cs = "Malaysia", main = "Malaysia (with time settings)",
         var = "tariff", ylim = c(-10, 60))


###################################################
### code chunk number 59: mmap1
###################################################
missmap(a.out)


###################################################
### code chunk number 60: mmap2
###################################################
getOption("SweaveHooks")[["fig"]]()
missmap(a.out)


###################################################
### code chunk number 61: amelia.Rnw:1775-1777
###################################################
a.out <- transform(a.out, lgdp = log(gdp.pc))
head(a.out$imputations[[1]][,c("country", "year","gdp.pc", "lgdp")])


###################################################
### code chunk number 62: amelia.Rnw:1780-1781
###################################################
a.out <- transform(a.out, pol_gdp = polity * gdp.pc)


###################################################
### code chunk number 63: amelia.Rnw:1785-1786
###################################################
summary(a.out)


###################################################
### code chunk number 64: amelia.Rnw:1803-1806
###################################################
require(Zelig)
z5 <- zls$new()
z5$zelig(tariff ~ polity + pop + gdp.pc + year + country, data = freetrade)


###################################################
### code chunk number 65: amelia.Rnw:1809-1810
###################################################
z5


###################################################
### code chunk number 66: amelia.Rnw:1817-1819
###################################################
z5 <- zls$new()
z5$zelig(tariff ~ polity + pop + gdp.pc + year + country, data = a.out)


###################################################
### code chunk number 67: amelia.Rnw:1822-1823
###################################################
z5


###################################################
### code chunk number 68: amelia.Rnw:1843-1853
###################################################
b.out<-NULL
se.out<-NULL
for(i in 1:a.out$m) {
  ols.out <- lm(tariff ~ polity + pop + gdp.pc, data = a.out$imputations[[i]])
  b.out <- rbind(b.out, ols.out$coef)
  se.out <- rbind(se.out, coef(summary(ols.out))[,2])
}

combined.results <- mi.meld(q = b.out, se = se.out)
combined.results


