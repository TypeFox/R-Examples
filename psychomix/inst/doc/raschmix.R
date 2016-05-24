### R code from vignette source 'raschmix.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(width = 70, prompt = "R> ", continue = "+  ")
library("psychomix")
data("VerbalAggression", package = "psychotools")
set.seed(1090)
cache <- FALSE


###################################################
### code chunk number 2: data-rost2
###################################################
set.seed(1)
r2 <- simRaschmix(design = "rost2")
d <- data.frame(
  x1 = rbinom(nrow(r2), prob = c(0.4, 0.6)[attr(r2, "cluster")], size = 1),
  x2 = rnorm(nrow(r2))
)
d$resp <- r2


###################################################
### code chunk number 3: fitRost0 (eval = FALSE)
###################################################
## set.seed(2)
## m1 <- raschmix(r2, k = 1:3)
## m1


###################################################
### code chunk number 4: fitRost
###################################################
if(cache & file.exists("m1.rda")) {
  load("m1.rda")
} else {
set.seed(2)
m1 <- raschmix(r2, k = 1:3)
m1
if(cache) {
  save(m1, file = "m1.rda")
} else {
  if(file.exists("m1.rda")) file.remove("m1.rda")
}
}


###################################################
### code chunk number 5: printRost
###################################################
m1


###################################################
### code chunk number 6: plotRost
###################################################
plot(m1)


###################################################
### code chunk number 7: selectRost
###################################################
BIC(m1)
m1b <- getModel(m1, which = "BIC")
summary(m1b)


###################################################
### code chunk number 8: itemParaRost
###################################################
parameters(m1b, "item")
worth(m1b)
attr(r2, "difficulty")


###################################################
### code chunk number 9: clustersRost
###################################################
table(model = clusters(m1b), true = attr(r2, "cluster"))


###################################################
### code chunk number 10: fitMeanVar0 (eval = FALSE)
###################################################
## set.seed(3)
## m2 <- raschmix(data = r2, k = 1:3, scores = "meanvar")


###################################################
### code chunk number 11: fitMeanVar
###################################################
if(cache & file.exists("m2.rda")) {
  load("m2.rda")
} else {
set.seed(3)
m2 <- raschmix(data = r2, k = 1:3, scores = "meanvar")
if(cache) {
  save(m2, file = "m2.rda")
} else {
  if(file.exists("m2.rda")) file.remove("m2.rda")
}
}


###################################################
### code chunk number 12: selectMeanVar
###################################################
m2
m2b <- getModel(m2, which = "BIC")


###################################################
### code chunk number 13: itemCompPlot
###################################################
par(mfrow = c(1,2))
plot(m1b, pos = "top")
for(i in 1:2) lines(attr(r2, "difficulty")[,i], lty = 2, type = "b")
plot(m2b, pos = "top")
for(i in 1:2) lines(attr(r2, "difficulty")[,i], lty = 2, type = "b")


###################################################
### code chunk number 14: compNPara
###################################################
logLik(m2b)
logLik(m1b)


###################################################
### code chunk number 15: scoreProbs:Rost2
###################################################
parameters(m2b, which = "score")
scoreProbs(m2b)


###################################################
### code chunk number 16: fitMeanVarConc0 (eval = FALSE)
###################################################
## set.seed(4)
## cm2 <- raschmix(resp ~ x1 + x2, data = d, k = 1:3, scores = "meanvar")


###################################################
### code chunk number 17: fitMeanVarConc
###################################################
if(cache & file.exists("cm2.rda")) {
  load("cm2.rda")
} else {
set.seed(4)
cm2 <- raschmix(resp ~ x1 + x2, data = d, k = 1:3, scores = "meanvar")
if(cache) {
  save(cm2, file = "cm2.rda")
} else {
  if(file.exists("cm2.rda")) file.remove("cm2.rda")
}
}


###################################################
### code chunk number 18: selectMeanVarConc
###################################################
rbind(m2 = BIC(m2), cm2 = BIC(cm2))


###################################################
### code chunk number 19: pValueMeanVarConc
###################################################
cm2b <- getModel(cm2, which = "BIC")
tStat <- 2 * (logLik(cm2b) - logLik(m2b))
pValue <- pchisq(tStat, attr(logLik(cm2b), "df") - attr(logLik(m2b), "df"), lower.tail = FALSE)
if(pValue < 0.001) pValue <- "< 0.001"


###################################################
### code chunk number 20: concPara
###################################################
cm2b <- getModel(cm2, which = "BIC")
parameters(cm2b, which = "concomitant")


###################################################
### code chunk number 21: concTable
###################################################
table(x1 = d$x1, clusters = clusters(cm2b))


###################################################
### code chunk number 22: dataVerbal
###################################################
data("VerbalAggression", package = "psychotools")
VerbalAggression$resp2 <- VerbalAggression$resp2[, 1:12]
va12 <- subset(VerbalAggression,
  rowSums(resp2) > 0 & rowSums(resp2) < 12)
colnames(va12$resp2)


###################################################
### code chunk number 23: vMeanVarFit0 (eval = FALSE)
###################################################
## set.seed(1)
## va12_mix1 <- raschmix(resp2 ~ 1, data = va12, k = 1:4, scores = "meanvar")
## set.seed(2)
## va12_mix2 <- raschmix(resp2 ~ gender + anger, data = va12, k = 1:4,
##   scores = "meanvar")


###################################################
### code chunk number 24: vMeanVarFit
###################################################
if(cache & file.exists("va12_mix.rda")) {
  load("va12_mix.rda")
} else {
set.seed(1)
va12_mix1 <- raschmix(resp2 ~ 1, data = va12, k = 1:4, scores = "meanvar")
set.seed(2)
va12_mix2 <- raschmix(resp2 ~ gender + anger, data = va12, k = 1:4,
  scores = "meanvar")
if(cache) {
  save(va12_mix1, va12_mix2, file = "va12_mix.rda")
} else {
  if(file.exists("va12_mix.rda")) file.remove("va12_mix.rda")
}
}


###################################################
### code chunk number 25: selectVMeanVar
###################################################
rbind(BIC(va12_mix1), BIC(va12_mix2))
va12_mix3 <- getModel(va12_mix2, which = "3")


###################################################
### code chunk number 26: VMeanVarPValue
###################################################
va12_mix1b <- getModel(va12_mix1, which = "3")
va12_mix2b <- getModel(va12_mix2, which = "3")

tStatVA <- 2 * (logLik(va12_mix2b) - logLik(va12_mix1b))
pValueVA <- pchisq(tStatVA, attr(logLik(va12_mix2b), "df") - 
   attr(logLik(va12_mix1b), "df"), lower.tail = FALSE)
if(pValueVA < 0.001) pValueVA <- "< 0.001"


###################################################
### code chunk number 27: plotVroot
###################################################
trellis.par.set(theme = standard.theme(color = FALSE))
print(histogram(va12_mix3))


###################################################
### code chunk number 28: plotVMeanVar
###################################################
trellis.par.set(theme = standard.theme(color = FALSE))
print(xyplot(va12_mix3))


###################################################
### code chunk number 29: concomitantsV
###################################################
parameters(getModel(va12_mix2, which = "3"), which = "concomitant")


###################################################
### code chunk number 30: appendixFit0 (eval = FALSE)
###################################################
## set.seed(4)
## fcm2 <- stepFlexmix(resp ~ 1, data = d, k = 1:3,
##   model = FLXMCrasch(scores = "meanvar"),
##   concomitant = FLXPmultinom(~ x1 + x2))


###################################################
### code chunk number 31: appendixFit
###################################################
if(cache & file.exists("fcm2.rda")) {
  load("fcm2.rda")
} else {
set.seed(4)
fcm2 <- stepFlexmix(resp ~ 1, data = d, k = 1:3,
  model = FLXMCrasch(scores = "meanvar"),
  concomitant = FLXPmultinom(~ x1 + x2))
if(cache) {
  save(fcm2, file = "fcm2.rda")
} else {
  if(file.exists("fcm2.rda")) file.remove("fcm2.rda")
}
}


###################################################
### code chunk number 32: appendixSelectModel
###################################################
rbind(cm2 = BIC(cm2), fcm2 = BIC(fcm2))
fcm2b <- getModel(fcm2, which = "BIC")


###################################################
### code chunk number 33: appendixParametersConc
###################################################
cbind(parameters(cm2b, which = "concomitant"),
  parameters(fcm2b, which = "concomitant"))


###################################################
### code chunk number 34: appendixParametes
###################################################
parameters(fcm2b, which = "model")


