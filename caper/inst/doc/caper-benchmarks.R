### R code from vignette source 'caper-benchmarks.rnw'

###################################################
### code chunk number 1: caper-benchmarks.rnw:23-25
###################################################
options(width=85)
library(xtable)


###################################################
### code chunk number 2: Datasets
###################################################
library(caper)
library(xtable)
data(benchTestInputs)
test <- comparative.data(testTree, testData, tips, vcv=TRUE, vcv.dim=3)
benchDicho <- comparative.data(benchTreeDicho, benchData, node, na.omit=FALSE)
benchPoly <- comparative.data(benchTreePoly, benchData, node, na.omit=FALSE)
print(test)
print(benchDicho)
print(benchPoly)


###################################################
### code chunk number 3: pglsTests (eval = FALSE)
###################################################
## bnds <- list(lambda=c(0.01,1), kappa=c(0.01,3), delta=c(0.01,3))
## 
## bnds <- list(lambda=c(0,1), kappa=c(0,3), delta=c(0,3))
## nul <- pgls(V1 ~ V2, test, bounds=bnds, lambda=0, kappa=0, delta=0)
## fix <- pgls(V1 ~ V2, test, bounds=bnds, lambda=0.5, kappa=0.5, delta=0.5)
## 
## kld <- pgls(V1 ~ V2, test, lambda=1, delta=1, kappa=1)
## Kld <- pgls(V1 ~ V2, test, lambda=1, delta=1, kappa="ML")
## kLd <- pgls(V1 ~ V2, test, lambda="ML", delta=1, kappa=1)
## klD <- pgls(V1 ~ V2, test, lambda=1, delta="ML", kappa=1)
## kLD <- pgls(V1 ~ V2, test, lambda="ML", delta="ML", kappa=1)
## KlD <- pgls(V1 ~ V2, test, lambda=1, delta="ML", kappa="ML")
## KLd <- pgls(V1 ~ V2, test, lambda="ML", delta=1, kappa="ML")
## KLD <- pgls(V1 ~ V2, test, lambda="ML", delta="ML", kappa="ML")
## 


###################################################
### code chunk number 4: pglsCheatyLoad
###################################################
data(pglsBenchmarks)


###################################################
### code chunk number 5: pglsMerge
###################################################
pglsMods <- list(nul=nul, fix=fix, kld = kld, Kld = Kld, kLd = kLd, 
	             klD = klD, kLD = kLD, KlD = KlD, KLd = KLd, KLD = KLD)

pglslogLik <- sapply(pglsMods, logLik)
pglsCoefs <- t(sapply(pglsMods, coef))
pglsSigma <- sapply(pglsMods, '[[', 'RMS')
pglsR.sq <- sapply(pglsMods, function(X) summary(X)$r.squared)
pglsParam <- t(sapply(pglsMods, '[[', 'param'))

pglsModTab <-  data.frame(mods=names(pglsMods), pglslogLik, pglsCoefs, pglsSigma, pglsR.sq, pglsParam)

data(benchBayesTraitsOutputs)


###################################################
### code chunk number 6: pglsTabRes
###################################################
pglsXt <- xtable(pglsModTab)
print(pglsXt, only.contents=TRUE, include.colnames=FALSE, include.rownames=FALSE)


###################################################
### code chunk number 7: bayesTabRes
###################################################
bayesXt <- xtable(BayesTraitsMods)
print(bayesXt, only.contents=TRUE, include.colnames=FALSE, include.rownames=FALSE)


###################################################
### code chunk number 8: crunchRun
###################################################
crunch.CrDi657 <- crunch(contRespNA ~ contExp1NA + contExp2NA,  data=benchDicho, polytomy.brlen=1)
crunch.CrDi213 <- crunch(contResp ~ contExp1 + contExp2, data=benchDicho, polytomy.brlen=1)
crunch.CrPl413 <- crunch(contResp ~ contExp1NoVar + contExp2, data=benchPoly, polytomy.brlen=1)
crunch.CrPl213 <- crunch(contResp ~ contExp1 + contExp2, data=benchPoly, polytomy.brlen=1)
crunch.CrPl657 <- crunch(contRespNA ~ contExp1NA + contExp2NA, data=benchPoly, polytomy.brlen=1)


###################################################
### code chunk number 9: crunchMerge
###################################################
crunch.CrDi657.tab <- caic.table(crunch.CrDi657, CAIC.codes=TRUE)
crunch.CrDi213.tab <- caic.table(crunch.CrDi213, CAIC.codes=TRUE)
crunch.CrPl413.tab <- caic.table(crunch.CrPl413, CAIC.codes=TRUE)
crunch.CrPl213.tab <- caic.table(crunch.CrPl213, CAIC.codes=TRUE)
crunch.CrPl657.tab <- caic.table(crunch.CrPl657, CAIC.codes=TRUE)

data(benchCrunchOutputs)

crunch.CrDi213.tab <- merge(crunch.CrDi213.tab,  CAIC.CrDi213, by.x="CAIC.code", by.y="Code", suffixes=c(".crunch", ".CAIC"))
crunch.CrDi657.tab <- merge(crunch.CrDi657.tab,  CAIC.CrDi657, by.x="CAIC.code", by.y="Code", suffixes=c(".crunch", ".CAIC"))
crunch.CrPl213.tab <- merge(crunch.CrPl213.tab,  CAIC.CrPl213, by.x="CAIC.code", by.y="Code", suffixes=c(".crunch", ".CAIC"))
crunch.CrPl413.tab <- merge(crunch.CrPl413.tab,  CAIC.CrPl413, by.x="CAIC.code", by.y="Code", suffixes=c(".crunch", ".CAIC"))
crunch.CrPl657.tab <- merge(crunch.CrPl657.tab,  CAIC.CrPl657, by.x="CAIC.code", by.y="Code", suffixes=c(".crunch", ".CAIC"))


###################################################
### code chunk number 10: crunchPlot
###################################################


par(mfrow=c(5,2), mar=c(2,2,1,1), mgp=c(1,0,0), tcl=-.1)

plot(contResp.CAIC ~ contExp1.CAIC, data=crunch.CrDi213.tab, xlab="contExp1", ylab="contResp")
with(crunch.CrDi213.tab, points(x=contExp1.crunch, y=contResp.crunch, col="red", pch=3, cex=0.8))

plot(contResp.CAIC ~ contExp2.CAIC, data=crunch.CrDi213.tab, xlab="contExp2", ylab="contResp")
with(crunch.CrDi213.tab, points(x=contExp2.crunch, y=contResp.crunch, col="red", pch=3, cex=0.8))


plot(contRespNA.CAIC ~ contExp1NA.CAIC, data=crunch.CrDi657.tab, xlab="contExp1NA", ylab="contRespNA")
with(crunch.CrDi657.tab, points(x=contExp1NA.crunch, y=contRespNA.crunch, col="red", pch=3, cex=0.8))

plot(contRespNA.CAIC ~ contExp2NA.CAIC, data=crunch.CrDi657.tab, xlab="contExp2NA", ylab="contResp")
with(crunch.CrDi657.tab, points(x=contExp2NA.crunch, y=contRespNA.crunch, col="red", pch=3, cex=0.8))


plot(contResp.CAIC ~ contExp1.CAIC, data=crunch.CrPl213.tab, xlab="contExp1", ylab="contResp")
with(crunch.CrPl213.tab, points(x=contExp1.crunch, y=contResp.crunch, col="red", pch=3, cex=0.8))

plot(contResp.CAIC ~ contExp2.CAIC, data=crunch.CrPl213.tab, xlab="contExp2", ylab="contResp")
with(crunch.CrPl213.tab, points(x=contExp2.crunch, y=contResp.crunch, col="red", pch=3, cex=0.8))


plot(contResp.CAIC ~ contExp1NoVar.CAIC, data=crunch.CrPl413.tab, xlab="contExp1NoVar", ylab="contResp")
with(crunch.CrPl413.tab, points(x=contExp1NoVar.crunch, y=contResp.crunch, col="red", pch=3, cex=0.8))

plot(contResp.CAIC ~ contExp2.CAIC, data=crunch.CrPl413.tab, xlab="contExp2", ylab="contResp")
with(crunch.CrPl413.tab, points(x=contExp2.crunch, y=contResp.crunch, col="red", pch=3, cex=0.8))


plot(contRespNA.CAIC ~ contExp1NA.CAIC, data=crunch.CrPl657.tab, xlab="contExp1NA", ylab="contRespNA")
with(crunch.CrPl657.tab, points(x=contExp1NA.crunch, y=contRespNA.crunch, col="red", pch=3, cex=0.8))

plot(contRespNA.CAIC ~ contExp2NA.CAIC, data=crunch.CrPl657.tab, xlab="contExp2NA", ylab="contResp")
with(crunch.CrPl657.tab, points(x=contExp2NA.crunch, y=contRespNA.crunch, col="red", pch=3, cex=0.8))




###################################################
### code chunk number 11: crunchCor
###################################################

# create a table of the correlations 
diffRangeTab <- data.frame(analysis=c("CrDi213","CrDi657","CrPl213","CrPl413","CrPl657"),
						RespCor=numeric(5), Exp1Cor=numeric(5), Exp2Cor=numeric(5))

diffRangeTab[1, 2:4] <- diag(cor(crunch.CrDi213.tab[,2:4], crunch.CrDi213.tab[,c(11,10,12)]))
diffRangeTab[2, 2:4] <- diag(cor(crunch.CrDi657.tab[,2:4], crunch.CrDi657.tab[,c(11,10,12)]))
diffRangeTab[3, 2:4] <- diag(cor(crunch.CrPl213.tab[,2:4], crunch.CrPl213.tab[,c(11,10,12)]))
diffRangeTab[4, 2:4] <- diag(cor(crunch.CrPl413.tab[,2:4], crunch.CrPl413.tab[,c(11,10,12)]))
diffRangeTab[5, 2:4] <- diag(cor(crunch.CrPl657.tab[,2:4], crunch.CrPl657.tab[,c(11,10,12)]))

diffRangeTab <- xtable(diffRangeTab, digits=4)
caption(diffRangeTab) <- "Correlations between values of crunch and CAIC contrasts."
label(diffRangeTab) <- "crunchCorr"
print(diffRangeTab, include.rownames=FALSE)


###################################################
### code chunk number 12: brunchRun
###################################################
brunch.BrDi813  <- brunch(contResp ~ binFact + contExp2, data=benchDicho)
brunch.BrDi913  <- brunch(contResp ~ triFact + contExp2, data=benchDicho)
brunch.BrDi1057 <- brunch(contRespNA ~ binFactNA + contExp2NA, data=benchDicho)
brunch.BrDi1157 <- brunch(contRespNA ~ triFactNA + contExp2NA, data=benchDicho)
brunch.BrPl813  <- brunch(contResp ~ binFact + contExp2, data=benchPoly)
brunch.BrPl913  <- brunch(contResp ~ triFact + contExp2, data=benchPoly)
brunch.BrPl1057 <- brunch(contRespNA ~ binFactNA + contExp2NA, data=benchPoly)
brunch.BrPl1157 <- brunch(contRespNA ~ triFactNA + contExp2NA, data=benchPoly)



###################################################
### code chunk number 13: brunchMerge
###################################################
brunch.BrDi813.tab  <- caic.table(brunch.BrDi813, CAIC.codes=TRUE)
brunch.BrDi913.tab  <- caic.table(brunch.BrDi913, CAIC.codes=TRUE)
brunch.BrDi1057.tab <- caic.table(brunch.BrDi1057, CAIC.codes=TRUE)
brunch.BrDi1157.tab <- caic.table(brunch.BrDi1157, CAIC.codes=TRUE)
brunch.BrPl813.tab  <- caic.table(brunch.BrPl813, CAIC.codes=TRUE)
brunch.BrPl913.tab  <- caic.table(brunch.BrPl913, CAIC.codes=TRUE)
brunch.BrPl1057.tab <- caic.table(brunch.BrPl1057, CAIC.codes=TRUE)
brunch.BrPl1157.tab <- caic.table(brunch.BrPl1157, CAIC.codes=TRUE)

data(benchBrunchOutputs)

brunch.BrDi813.tab  <- merge(brunch.BrDi813.tab , CAIC.BrDi813 , by.x="CAIC.code", by.y="Code", suffixes=c(".brunch", ".CAIC"))
brunch.BrDi913.tab  <- merge(brunch.BrDi913.tab , CAIC.BrDi913 , by.x="CAIC.code", by.y="Code", suffixes=c(".brunch", ".CAIC"))
brunch.BrDi1057.tab <- merge(brunch.BrDi1057.tab, CAIC.BrDi1057, by.x="CAIC.code", by.y="Code", suffixes=c(".brunch", ".CAIC"))
brunch.BrDi1157.tab <- merge(brunch.BrDi1157.tab, CAIC.BrDi1157, by.x="CAIC.code", by.y="Code", suffixes=c(".brunch", ".CAIC"))
brunch.BrPl813.tab  <- merge(brunch.BrPl813.tab , CAIC.BrPl813 , by.x="CAIC.code", by.y="Code", suffixes=c(".brunch", ".CAIC"))
brunch.BrPl913.tab  <- merge(brunch.BrPl913.tab , CAIC.BrPl913 , by.x="CAIC.code", by.y="Code", suffixes=c(".brunch", ".CAIC"))
brunch.BrPl1057.tab <- merge(brunch.BrPl1057.tab, CAIC.BrPl1057, by.x="CAIC.code", by.y="Code", suffixes=c(".brunch", ".CAIC"))
brunch.BrPl1157.tab <- merge(brunch.BrPl1157.tab, CAIC.BrPl1157, by.x="CAIC.code", by.y="Code", suffixes=c(".brunch", ".CAIC"))


###################################################
### code chunk number 14: brunchBin
###################################################
par(mfrow=c(4,2), mar=c(2,2,1,1), mgp=c(1,0,0), tcl=-.1, oma=c(0,3,0,0))

plot(contResp.CAIC ~ binFact.CAIC, data=brunch.BrDi813.tab)
points(contResp.brunch ~ binFact.brunch, data=brunch.BrDi813.tab, col="red", pch=3, cex=0.8)
mtext('BrDi813', side=2, line=3)

plot(contResp.CAIC ~ contExp2.CAIC, data=brunch.BrDi813.tab)
points(contResp.brunch ~ contExp2.brunch, data=brunch.BrDi813.tab,col="red", pch=3, cex=0.8)


plot(contRespNA.CAIC ~ binFactNA.CAIC, data=brunch.BrDi1057.tab)
points(contRespNA.brunch ~ binFactNA.brunch, data=brunch.BrDi1057.tab, col="red", pch=3, cex=0.8)
mtext('BrDi1057', side=2, line=3)

plot(contRespNA.CAIC ~ contExp2NA.CAIC, data=brunch.BrDi1057.tab)
points(contRespNA.brunch ~ contExp2NA.brunch, data=brunch.BrDi1057.tab,col="red", pch=3, cex=0.8)

plot(contResp.CAIC ~ binFact.CAIC, data=brunch.BrPl813.tab)
points(contResp.brunch ~ binFact.brunch, data=brunch.BrPl813.tab, col="red", pch=3, cex=0.8)
mtext('BrPl813', side=2, line=3)

plot(contResp.CAIC ~ contExp2.CAIC, data=brunch.BrPl813.tab)
points(contResp.brunch ~ contExp2.brunch, data=brunch.BrPl813.tab,col="red", pch=3, cex=0.8)


plot(contRespNA.CAIC ~ binFactNA.CAIC, data=brunch.BrPl1057.tab)
points(contRespNA.brunch ~ binFactNA.brunch, data=brunch.BrPl1057.tab, col="red", pch=3, cex=0.8)
mtext('BrPl1057', side=2, line=3)

plot(contRespNA.CAIC ~ contExp2NA.CAIC, data=brunch.BrPl1057.tab)
points(contRespNA.brunch ~ contExp2NA.brunch, data=brunch.BrPl1057.tab,col="red", pch=3, cex=0.8)


###################################################
### code chunk number 15: brunchTri
###################################################
par(mfrow=c(4,2), mar=c(2,2,1,1), mgp=c(1,0,0), tcl=-.1, oma=c(0,3,0,0))

plot(contResp.CAIC ~ triFact.CAIC, data=brunch.BrDi913.tab)
points(contResp.brunch ~ triFact.brunch, data=brunch.BrDi913.tab, col="red", pch=3, cex=0.8)
mtext('BrDi913', side=2, line=3)

plot(contResp.CAIC ~ contExp2.CAIC, data=brunch.BrDi913.tab)
points(contResp.brunch ~ contExp2.brunch, data=brunch.BrDi913.tab,col="red", pch=3, cex=0.8)


plot(contRespNA.CAIC ~triFactNA.CAIC, data=brunch.BrDi1157.tab)
points(contRespNA.brunch ~ triFactNA.brunch, data=brunch.BrDi1157.tab, col="red", pch=3, cex=0.8)
mtext('BrDi1157', side=2, line=3)

plot(contRespNA.CAIC ~ contExp2NA.CAIC, data=brunch.BrDi1157.tab)
points(contRespNA.brunch ~ contExp2NA.brunch, data=brunch.BrDi1157.tab,col="red", pch=3, cex=0.8)

plot(contResp.CAIC ~ triFact.CAIC, data=brunch.BrPl913.tab)
points(contResp.brunch ~ triFact.brunch, data=brunch.BrPl913.tab, col="red", pch=3, cex=0.8)
mtext('BrPl913', side=2, line=3)

plot(contResp.CAIC ~ contExp2.CAIC, data=brunch.BrPl913.tab)
points(contResp.brunch ~ contExp2.brunch, data=brunch.BrPl913.tab,col="red", pch=3, cex=0.8)


plot(contRespNA.CAIC ~triFactNA.CAIC, data=brunch.BrPl1157.tab)
points(contRespNA.brunch ~ triFactNA.brunch, data=brunch.BrPl1157.tab, col="red", pch=3, cex=0.8)
mtext('BrPl1157', side=2, line=3)

plot(contRespNA.CAIC ~ contExp2NA.CAIC, data=brunch.BrPl1157.tab)
points(contRespNA.brunch ~ contExp2NA.brunch, data=brunch.BrPl1157.tab,col="red", pch=3, cex=0.8)




###################################################
### code chunk number 16: brunchCor
###################################################
# create a table of the correlations

diffRangeTab <- data.frame(analysis=c("BrDi813","BrPl813","BrDi1057","BrPl1057",
                                      "BrDi913","BrPl913","BrDi1157","BrPl1157"),
						   RespCor=numeric(8), Exp1Cor=numeric(8), Exp2Cor=numeric(8))

diffRangeTab[1, 2:4] <- diag(cor(brunch.BrDi813.tab[,2:4], brunch.BrDi813.tab[,c(11,10,12)]))
diffRangeTab[2, 2:4] <- diag(cor(brunch.BrDi1057.tab[,2:4], brunch.BrDi1057.tab[,c(11,10,12)]))
diffRangeTab[3, 2:4] <- diag(cor(brunch.BrDi913.tab[,2:4], brunch.BrDi913.tab[,c(11,10,12)]))
diffRangeTab[4, 2:4] <- diag(cor(brunch.BrDi1157.tab[,2:4], brunch.BrDi1157.tab[,c(11,10,12)]))
diffRangeTab[5, 2:4] <- diag(cor(brunch.BrPl813.tab[,2:4], brunch.BrPl813.tab[,c(11,10,12)]))
diffRangeTab[6, 2:4] <- diag(cor(brunch.BrPl1057.tab[,2:4], brunch.BrPl1057.tab[,c(11,10,12)]))
diffRangeTab[7, 2:4] <- diag(cor(brunch.BrPl913.tab[,2:4], brunch.BrPl913.tab[,c(11,10,12)]))
diffRangeTab[8, 2:4] <- diag(cor(brunch.BrPl1157.tab[,2:4], brunch.BrPl1157.tab[,c(11,10,12)]))

diffRangeTab <- xtable(diffRangeTab, digits=5)
caption(diffRangeTab) <- "Range in the differences between brunch and CAIC contrasts."
label(diffRangeTab) <- "brunchCorr"
print(diffRangeTab, include.rownames=FALSE)


###################################################
### code chunk number 17: macroRun
###################################################

DichSpp23.RRD <- macrocaic(sppRichTips ~ contExp1  + contExp2, macroMethod = "RRD", data=benchDicho)
PolySpp23.RRD <- macrocaic(sppRichTips ~ contExp1  + contExp2, macroMethod = "RRD", data=benchPoly)
DichTax23.RRD <- macrocaic(sppRichTaxa ~ contExp1  + contExp2, macroMethod = "RRD", data=benchDicho)
PolyTax23.RRD <- macrocaic(sppRichTaxa ~ contExp1  + contExp2, macroMethod = "RRD", data=benchPoly)
DichSpp67.RRD <- macrocaic(sppRichTips ~ contExp1NA + contExp2NA, macroMethod = "RRD", data=benchDicho)
PolySpp67.RRD <- macrocaic(sppRichTips ~ contExp1NA + contExp2NA, macroMethod = "RRD", data=benchPoly)
DichTax67.RRD <- macrocaic(sppRichTaxa ~ contExp1NA + contExp2NA, macroMethod = "RRD", data=benchDicho)
PolyTax67.RRD <- macrocaic(sppRichTaxa ~ contExp1NA + contExp2NA, macroMethod = "RRD", data=benchPoly)

DichSpp23.PDI <- macrocaic(sppRichTips ~ contExp1  + contExp2, macroMethod = "PDI", data=benchDicho)
PolySpp23.PDI <- macrocaic(sppRichTips ~ contExp1  + contExp2, macroMethod = "PDI", data=benchPoly)
DichTax23.PDI <- macrocaic(sppRichTaxa ~ contExp1  + contExp2, macroMethod = "PDI", data=benchDicho)
PolyTax23.PDI <- macrocaic(sppRichTaxa ~ contExp1  + contExp2, macroMethod = "PDI", data=benchPoly)
DichSpp67.PDI <- macrocaic(sppRichTips ~ contExp1NA + contExp2NA, macroMethod = "PDI", data=benchDicho)
PolySpp67.PDI <- macrocaic(sppRichTips ~ contExp1NA + contExp2NA, macroMethod = "PDI", data=benchPoly)
DichTax67.PDI <- macrocaic(sppRichTaxa ~ contExp1NA + contExp2NA, macroMethod = "PDI", data=benchDicho)
PolyTax67.PDI <- macrocaic(sppRichTaxa ~ contExp1NA + contExp2NA, macroMethod = "PDI", data=benchPoly)



###################################################
### code chunk number 18: macroMerge
###################################################
macro.DichSpp23 <- merge(caic.table(DichSpp23.RRD, CAIC.codes=TRUE), caic.table(DichSpp23.PDI, CAIC.codes=TRUE), by="CAIC.code", suffixes=c(".mcRRD", ".mcPDI"))
macro.PolySpp23 <- merge(caic.table(PolySpp23.RRD, CAIC.codes=TRUE), caic.table(PolySpp23.PDI, CAIC.codes=TRUE), by="CAIC.code", suffixes=c(".mcRRD", ".mcPDI"))
macro.DichTax23 <- merge(caic.table(DichTax23.RRD, CAIC.codes=TRUE), caic.table(DichTax23.PDI, CAIC.codes=TRUE), by="CAIC.code", suffixes=c(".mcRRD", ".mcPDI"))
macro.PolyTax23 <- merge(caic.table(PolyTax23.RRD, CAIC.codes=TRUE), caic.table(PolyTax23.PDI, CAIC.codes=TRUE), by="CAIC.code", suffixes=c(".mcRRD", ".mcPDI"))
macro.DichSpp67 <- merge(caic.table(DichSpp67.RRD, CAIC.codes=TRUE), caic.table(DichSpp67.PDI, CAIC.codes=TRUE), by="CAIC.code", suffixes=c(".mcRRD", ".mcPDI"))
macro.PolySpp67 <- merge(caic.table(PolySpp67.RRD, CAIC.codes=TRUE), caic.table(PolySpp67.PDI, CAIC.codes=TRUE), by="CAIC.code", suffixes=c(".mcRRD", ".mcPDI"))
macro.DichTax67 <- merge(caic.table(DichTax67.RRD, CAIC.codes=TRUE), caic.table(DichTax67.PDI, CAIC.codes=TRUE), by="CAIC.code", suffixes=c(".mcRRD", ".mcPDI"))
macro.PolyTax67 <- merge(caic.table(PolyTax67.RRD, CAIC.codes=TRUE), caic.table(PolyTax67.PDI, CAIC.codes=TRUE), by="CAIC.code", suffixes=c(".mcRRD", ".mcPDI"))

data(benchMacroCAICOutputs)

macro.DichSpp23 <- merge(macro.DichSpp23, MacroCAIC.DiSpp23, by.x="CAIC.code", by.y="Code")
macro.PolySpp23 <- merge(macro.PolySpp23, MacroCAIC.PolySpp23, by.x="CAIC.code", by.y="Code")
macro.DichTax23 <- merge(macro.DichTax23, MacroCAIC.DiTax23, by.x="CAIC.code", by.y="Code")
macro.PolyTax23 <- merge(macro.PolyTax23, MacroCAIC.PolyTax23, by.x="CAIC.code", by.y="Code")
macro.DichSpp67 <- merge(macro.DichSpp67, MacroCAIC.DiSpp67, by.x="CAIC.code", by.y="Code")
macro.PolySpp67 <- merge(macro.PolySpp67, MacroCAIC.PolySpp67, by.x="CAIC.code", by.y="Code")
macro.DichTax67 <- merge(macro.DichTax67, MacroCAIC.DiTax67, by.x="CAIC.code", by.y="Code")
macro.PolyTax67 <- merge(macro.PolyTax67, MacroCAIC.PolyTax67, by.x="CAIC.code", by.y="Code")



###################################################
### code chunk number 19: macroCorr
###################################################
# create a table to keep track of the differences 
MacroCAICTab <- data.frame(analysis=c("DichSpp23", "PolySpp23", "DichTax23", "PolyTax23", 
                                      "DichSpp67", "PolySpp67", "DichTax67", "PolyTax67"),
						   RichCor.RRD=integer(8), Exp1Cor.RRD=integer(8), Exp2Cor.RRD=integer(8),
						   RichCor.PDI=integer(8), Exp1Cor.PDI=integer(8), Exp2Cor.PDI=integer(8))

MacroCAICTab[1, 2:4] <- diag(cor(macro.DichSpp23[, 2:4], macro.DichSpp23[, c(23,16,17)]))
MacroCAICTab[2, 2:4] <- diag(cor(macro.PolySpp23[, 2:4], macro.PolySpp23[, c(23,16,17)]))
MacroCAICTab[3, 2:4] <- diag(cor(macro.DichTax23[, 2:4], macro.DichTax23[, c(23,16,17)]))
MacroCAICTab[4, 2:4] <- diag(cor(macro.PolyTax23[, 2:4], macro.PolyTax23[, c(23,16,17)]))
MacroCAICTab[5, 2:4] <- diag(cor(macro.DichSpp67[, 2:4], macro.DichSpp67[, c(23,16,17)]))
MacroCAICTab[6, 2:4] <- diag(cor(macro.PolySpp67[, 2:4], macro.PolySpp67[, c(23,16,17)]))
MacroCAICTab[7, 2:4] <- diag(cor(macro.DichTax67[, 2:4], macro.DichTax67[, c(23,16,17)]))
MacroCAICTab[8, 2:4] <- diag(cor(macro.PolyTax67[, 2:4], macro.PolyTax67[, c(23,16,17)]))

MacroCAICTab[1, 5:7] <- diag(cor(macro.DichSpp23[, 9:11], macro.DichSpp23[, c(22,16,17)]))
MacroCAICTab[2, 5:7] <- diag(cor(macro.PolySpp23[, 9:11], macro.PolySpp23[, c(22,16,17)]))
MacroCAICTab[3, 5:7] <- diag(cor(macro.DichTax23[, 9:11], macro.DichTax23[, c(22,16,17)]))
MacroCAICTab[4, 5:7] <- diag(cor(macro.PolyTax23[, 9:11], macro.PolyTax23[, c(22,16,17)]))
MacroCAICTab[5, 5:7] <- diag(cor(macro.DichSpp67[, 9:11], macro.DichSpp67[, c(22,16,17)]))
MacroCAICTab[6, 5:7] <- diag(cor(macro.PolySpp67[, 9:11], macro.PolySpp67[, c(22,16,17)]))
MacroCAICTab[7, 5:7] <- diag(cor(macro.DichTax67[, 9:11], macro.DichTax67[, c(22,16,17)]))
MacroCAICTab[8, 5:7] <- diag(cor(macro.PolyTax67[, 9:11], macro.PolyTax67[, c(22,16,17)]))

MacroCAICTab <- xtable(MacroCAICTab, digits=5)
caption(MacroCAICTab) <- "Correlations between brunch and CAIC contrasts."
label(MacroCAICTab) <- "macroCorr"
print(MacroCAICTab, include.rownames=FALSE)


###################################################
### code chunk number 20: RRDComp
###################################################
par(mfrow=c(4,2), mar=c(2,2,1,1), mgp=c(1,0,0), tcl=-.1, oma=c(0,3,0,0))

plot(RRD ~ contExp1, data=macro.DichSpp23)
points(sppRichTips.mcRRD ~ contExp1.mcRRD, data=macro.DichSpp23, col="red", pch=3, cex=0.8)
mtext("DichSpp23", side=2, line=3)

plot(RRD ~ contExp2, data=macro.DichSpp23)
points(sppRichTips.mcRRD ~ contExp2.mcRRD, data=macro.DichSpp23, col="red", pch=3, cex=0.8)

plot(RRD ~ contExp1, data=macro.PolySpp23)
points(sppRichTips.mcRRD ~ contExp1.mcRRD, data=macro.PolySpp23, col="red", pch=3, cex=0.8)
mtext("PolySpp23", side=2, line=3)

plot(RRD ~ contExp2, data=macro.PolySpp23)
points(sppRichTips.mcRRD ~ contExp2.mcRRD, data=macro.PolySpp23, col="red", pch=3, cex=0.8)

plot(RRD ~ contExp1, data=macro.DichTax23)
points(sppRichTaxa.mcRRD ~ contExp1.mcRRD, data=macro.DichTax23, col="red", pch=3, cex=0.8)
mtext("DichTax23", side=2, line=3)

plot(RRD ~ contExp2, data=macro.DichTax23)
points(sppRichTaxa.mcRRD ~ contExp2.mcRRD, data=macro.DichTax23, col="red", pch=3, cex=0.8)

plot(RRD ~ contExp1, data=macro.PolyTax23)
points(sppRichTaxa.mcRRD ~ contExp1.mcRRD, data=macro.PolyTax23, col="red", pch=3, cex=0.8)
mtext("PolyTax23", side=2, line=3)

plot(RRD ~ contExp2, data=macro.PolyTax23)
points(sppRichTaxa.mcRRD ~ contExp2.mcRRD, data=macro.PolyTax23, col="red", pch=3, cex=0.8)


###################################################
### code chunk number 21: RRDMiss
###################################################
par(mfrow=c(4,2), mar=c(2,2,1,1), mgp=c(1,0,0), tcl=-.1, oma=c(0,3,0,0))

plot(RRD ~ contExp1NA, data=macro.DichSpp67)
points(sppRichTips.mcRRD ~ contExp1NA.mcRRD, data=macro.DichSpp67, col="red", pch=3, cex=0.8)
mtext("DichSpp67", side=2, line=3)

plot(RRD ~ contExp2NA, data=macro.DichSpp67)
points(sppRichTips.mcRRD ~ contExp2NA.mcRRD, data=macro.DichSpp67, col="red", pch=3, cex=0.8)

plot(RRD ~ contExp1NA, data=macro.PolySpp67)
points(sppRichTips.mcRRD ~ contExp1NA.mcRRD, data=macro.PolySpp67, col="red", pch=3, cex=0.8)
mtext("PolySpp67", side=2, line=3)

plot(RRD ~ contExp2NA, data=macro.PolySpp67)
points(sppRichTips.mcRRD ~ contExp2NA.mcRRD, data=macro.PolySpp67, col="red", pch=3, cex=0.8)

plot(RRD ~ contExp1NA, data=macro.DichTax67)
points(sppRichTaxa.mcRRD ~ contExp1NA.mcRRD, data=macro.DichTax67, col="red", pch=3, cex=0.8)
mtext("DichTax67", side=2, line=3)

plot(RRD ~ contExp2NA, data=macro.DichTax67)
points(sppRichTaxa.mcRRD ~ contExp2NA.mcRRD, data=macro.DichTax67, col="red", pch=3, cex=0.8)

plot(RRD ~ contExp1NA, data=macro.PolyTax67)
points(sppRichTaxa.mcRRD ~ contExp1NA.mcRRD, data=macro.PolyTax67, col="red", pch=3, cex=0.8)
mtext("PolyTax67", side=2, line=3)

plot(RRD ~ contExp2NA, data=macro.PolyTax67)
points(sppRichTaxa.mcRRD ~ contExp2NA.mcRRD, data=macro.PolyTax67, col="red", pch=3, cex=0.8)


###################################################
### code chunk number 22: PDIComp
###################################################
par(mfrow=c(4,2), mar=c(2,2,1,1), mgp=c(1,0,0), tcl=-.1, oma=c(0,3,0,0))

plot(PDI ~ contExp1, data=macro.DichSpp23)
points(sppRichTips.mcPDI ~ contExp1.mcPDI, data=macro.DichSpp23, col="red", pch=3, cex=0.8)
mtext("DichSpp23", side=2, line=3)

plot(PDI ~ contExp2, data=macro.DichSpp23)
points(sppRichTips.mcPDI ~ contExp2.mcPDI, data=macro.DichSpp23, col="red", pch=3, cex=0.8)

plot(PDI ~ contExp1, data=macro.PolySpp23)
points(sppRichTips.mcPDI ~ contExp1.mcPDI, data=macro.PolySpp23, col="red", pch=3, cex=0.8)
mtext("PolySpp23", side=2, line=3)

plot(PDI ~ contExp2, data=macro.PolySpp23)
points(sppRichTips.mcPDI ~ contExp2.mcPDI, data=macro.PolySpp23, col="red", pch=3, cex=0.8)

plot(PDI ~ contExp1, data=macro.DichTax23)
points(sppRichTaxa.mcPDI ~ contExp1.mcPDI, data=macro.DichTax23, col="red", pch=3, cex=0.8)
mtext("DichTax23", side=2, line=3)

plot(PDI ~ contExp2, data=macro.DichTax23)
points(sppRichTaxa.mcPDI ~ contExp2.mcPDI, data=macro.DichTax23, col="red", pch=3, cex=0.8)

plot(PDI ~ contExp1, data=macro.PolyTax23)
points(sppRichTaxa.mcPDI ~ contExp1.mcPDI, data=macro.PolyTax23, col="red", pch=3, cex=0.8)
mtext("PolyTax23", side=2, line=3)

plot(PDI ~ contExp2, data=macro.PolyTax23)
points(sppRichTaxa.mcPDI ~ contExp2.mcPDI, data=macro.PolyTax23, col="red", pch=3, cex=0.8)


###################################################
### code chunk number 23: PDIMiss
###################################################
par(mfrow=c(4,2), mar=c(2,2,1,1), mgp=c(1,0,0), tcl=-.1, oma=c(0,3,0,0))

plot(PDI ~ contExp1NA, data=macro.DichSpp67)
points(sppRichTips.mcPDI ~ contExp1NA.mcPDI, data=macro.DichSpp67, col="red", pch=3, cex=0.8)
mtext("DichSpp67", side=2, line=3)

plot(PDI ~ contExp2NA, data=macro.DichSpp67)
points(sppRichTips.mcPDI ~ contExp2NA.mcPDI, data=macro.DichSpp67, col="red", pch=3, cex=0.8)

plot(PDI ~ contExp1NA, data=macro.PolySpp67)
points(sppRichTips.mcPDI ~ contExp1NA.mcPDI, data=macro.PolySpp67, col="red", pch=3, cex=0.8)
mtext("PolySpp67", side=2, line=3)

plot(PDI ~ contExp2NA, data=macro.PolySpp67)
points(sppRichTips.mcPDI ~ contExp2NA.mcPDI, data=macro.PolySpp67, col="red", pch=3, cex=0.8)

plot(PDI ~ contExp1NA, data=macro.DichTax67)
points(sppRichTaxa.mcPDI ~ contExp1NA.mcPDI, data=macro.DichTax67, col="red", pch=3, cex=0.8)
mtext("DichTax67", side=2, line=3)

plot(PDI ~ contExp2NA, data=macro.DichTax67)
points(sppRichTaxa.mcPDI ~ contExp2NA.mcPDI, data=macro.DichTax67, col="red", pch=3, cex=0.8)

plot(PDI ~ contExp1NA, data=macro.PolyTax67)
points(sppRichTaxa.mcPDI ~ contExp1NA.mcPDI, data=macro.PolyTax67, col="red", pch=3, cex=0.8)
mtext("PolyTax67", side=2, line=3)

plot(PDI ~ contExp2NA, data=macro.PolyTax67)
points(sppRichTaxa.mcPDI ~ contExp2NA.mcPDI, data=macro.PolyTax67, col="red", pch=3, cex=0.8)



###################################################
### code chunk number 24: fuscoRunMerge
###################################################
fstest.DiSpp <- fusco.test(benchDicho, rich=sppRichTips)
fstest.DiTax <- fusco.test(benchDicho, rich=sppRichTaxa)
fstest.PlSpp <- fusco.test(benchPoly, rich=sppRichTips)
fstest.PlTax <- fusco.test(benchPoly, rich=sppRichTaxa)

fstest.DiSpp.Hist <- plot(fstest.DiSpp, I.prime=FALSE, plot=FALSE)
fstest.DiTax.Hist <- plot(fstest.DiTax, I.prime=FALSE, plot=FALSE)
fstest.PlSpp.Hist <- plot(fstest.PlSpp, I.prime=FALSE, plot=FALSE)
fstest.PlTax.Hist <- plot(fstest.PlTax, I.prime=FALSE, plot=FALSE)

data(benchFuscoOutputs)

DiSpp <- cbind(fstest.DiSpp.Hist, FuscoDiSpp$distTab) 
DiTax <- cbind(fstest.DiTax.Hist, FuscoDiTax$distTab)
PlSpp <- cbind(fstest.PlSpp.Hist, FuscoPolySpp$distTab)
PlTax <- cbind(fstest.PlTax.Hist, FuscoPolyTax$distTab)



###################################################
### code chunk number 25: fuscoPlots
###################################################
par(mfrow=c(2,2), mar=c(2,3,1,1), mgp=c(1,0,0), tcl=-.1, cex.axis=0.7)

barplot(t(DiSpp[,c(4,8)]), beside=TRUE, xlab='Corrected imbalance', ylab='Density')
axis(side=1, at=seq(0,30, by=3), label=seq(0, 1, by=0.1))
mtext(side=2, 'DiSpp', line=2)
barplot(t(DiTax[,c(4,8)]), beside=TRUE, xlab='Corrected imbalance', ylab='Density')
axis(side=1, at=seq(0,30, by=3), label=seq(0, 1, by=0.1))
mtext(side=2, 'DiTax', line=2)
barplot(t(PlSpp[,c(4,8)]), beside=TRUE, xlab='Corrected imbalance', ylab='Density')
axis(side=1, at=seq(0,30, by=3), label=seq(0, 1, by=0.1))
mtext(side=2, 'PlSpp', line=2)
barplot(t(PlTax[,c(4,8)]), beside=TRUE, xlab='Corrected imbalance', ylab='Density')
axis(side=1, at=seq(0,30, by=3), label=seq(0, 1, by=0.1))
mtext(side=2, 'PlTax', line=2)


###################################################
### code chunk number 26: caper-benchmarks.rnw:740-747
###################################################
data(syrphidae)

fstest.Syrph <- fusco.test(syrphidaeTree, dat=syrphidaeRich, rich=nSpp, names=genus)
summary(fstest.Syrph)
data(benchMesaOutputs)
mean(MeSA.I$Iprime)
median(MeSA.I$I)


