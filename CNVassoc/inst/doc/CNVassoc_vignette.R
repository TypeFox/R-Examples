### R code from vignette source 'CNVassoc_vignette.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: CNVassoc_vignette.Rnw:62-63
###################################################
options(width=120)


###################################################
### code chunk number 2: CNVassoc_vignette.Rnw:66-67
###################################################
library(CNVassoc)


###################################################
### code chunk number 3: CNVassoc_vignette.Rnw:71-72
###################################################
library(xtable)


###################################################
### code chunk number 4: CNVassoc_vignette.Rnw:94-96
###################################################
data(dataMLPA)
head(dataMLPA)


###################################################
### code chunk number 5: CNVassoc_vignette.Rnw:113-119
###################################################
par(mfrow=c(2,2),mar=c(3,4,3,1))
hist(dataMLPA$Gene1,main="Gene 1 signal histogram",xlab="",ylab="frequency")
hist(dataMLPA$Gene2,main="Gene 2 signal histogram",xlab="",ylab="frequency")
par(xaxs="i")
plot(density(dataMLPA$Gene1),main="Gene 1 signal density function",xlab="",ylab="density")
plot(density(dataMLPA$Gene2),main="Gene 2 signal density function",xlab="",ylab="density")


###################################################
### code chunk number 6: CNVassoc_vignette.Rnw:128-129
###################################################
plotSignal(dataMLPA$Gene2,case.control=dataMLPA$casco)


###################################################
### code chunk number 7: CNVassoc_vignette.Rnw:134-135
###################################################
plotSignal(dataMLPA$Gene2,case.control=dataMLPA$casco)


###################################################
### code chunk number 8: CNVassoc_vignette.Rnw:143-144
###################################################
plotSignal(dataMLPA$Gene2, case.control = dataMLPA$quanti)


###################################################
### code chunk number 9: CNVassoc_vignette.Rnw:150-151
###################################################
plotSignal(dataMLPA$Gene2,case.control=dataMLPA$quanti)


###################################################
### code chunk number 10: CNVassoc_vignette.Rnw:188-189
###################################################
cutpoints <- c(0.08470221, 0.40485249)


###################################################
### code chunk number 11: CNVassoc_vignette.Rnw:191-192
###################################################
cutpoints


###################################################
### code chunk number 12: CNVassoc_vignette.Rnw:226-228
###################################################
CNV.1 <- cnv(x = dataMLPA$Gene1, threshold.0 = 0.06, num.class = 3, 
mix.method = "mixdist")


###################################################
### code chunk number 13: CNVassoc_vignette.Rnw:243-244
###################################################
CNV.2 <- cnv(x = dataMLPA$Gene2, threshold.0 = 0.01, mix.method = "mixdist")


###################################################
### code chunk number 14: CNVassoc_vignette.Rnw:266-267
###################################################
probs.2 <- getProbs(CNV.2)


###################################################
### code chunk number 15: CNVassoc_vignette.Rnw:276-277
###################################################
CNV.2probs <- cnv(probs.2)


###################################################
### code chunk number 16: CNVassoc_vignette.Rnw:287-288
###################################################
CNV.2th <- cnv(x = dataMLPA$Gene2, cutoffs = cutpoints)


###################################################
### code chunk number 17: CNVassoc_vignette.Rnw:301-302
###################################################
CNV.1


###################################################
### code chunk number 18: CNVassoc_vignette.Rnw:309-310
###################################################
CNV.2


###################################################
### code chunk number 19: CNVassoc_vignette.Rnw:317-318
###################################################
CNV.2probs


###################################################
### code chunk number 20: CNVassoc_vignette.Rnw:325-331
###################################################
pdf("fig2a.pdf")
plot(CNV.1, case.control = dataMLPA$casco, main = "Gene 1")
dev.off()
pdf("fig2b.pdf")
plot(CNV.2, case.control = dataMLPA$casco, main = "Gene 2")
dev.off()


###################################################
### code chunk number 21: CNVassoc_vignette.Rnw:364-365
###################################################
plot(CNV.2probs, case.control=dataMLPA$casco)


###################################################
### code chunk number 22: CNVassoc_vignette.Rnw:391-393
###################################################
CNVassoc::getQualityScore(CNV.1, type = "class")
CNVassoc::getQualityScore(CNV.2, type = "class")


###################################################
### code chunk number 23: CNVassoc_vignette.Rnw:400-402
###################################################
CNVassoc::getQualityScore(CNV.1, type = "CNVtools")
CNVassoc::getQualityScore(CNV.2, type = "CNVtools")


###################################################
### code chunk number 24: CNVassoc_vignette.Rnw:409-411
###################################################
CNVassoc::getQualityScore(CNV.1, type = "CANARY")
CNVassoc::getQualityScore(CNV.2, type = "CANARY")


###################################################
### code chunk number 25: CNVassoc_vignette.Rnw:448-450
###################################################
model1mul <- CNVassoc(casco ~ CNV.1, data = dataMLPA, model = "mul")
model2mul <- CNVassoc(casco ~ CNV.2, data = dataMLPA, model = "mul")


###################################################
### code chunk number 26: CNVassoc_vignette.Rnw:461-462
###################################################
model1mul


###################################################
### code chunk number 27: CNVassoc_vignette.Rnw:468-469
###################################################
model2mul


###################################################
### code chunk number 28: CNVassoc_vignette.Rnw:484-485
###################################################
summary(model1mul)


###################################################
### code chunk number 29: CNVassoc_vignette.Rnw:492-493
###################################################
summary(model2mul)


###################################################
### code chunk number 30: CNVassoc_vignette.Rnw:505-508
###################################################
mod <- CNVassoc(quanti ~ CNV.2 + cov, family = "gaussian", 
data = dataMLPA, model = 'add', emsteps = 10)
mod


###################################################
### code chunk number 31: CNVassoc_vignette.Rnw:518-519
###################################################
summary(mod)


###################################################
### code chunk number 32: CNVassoc_vignette.Rnw:536-538
###################################################
CNVtest(model1mul, type = "Wald")
CNVtest(model1mul, type = "LRT")


###################################################
### code chunk number 33: CNVassoc_vignette.Rnw:545-547
###################################################
CNVtest(model2mul, type = "Wald")
CNVtest(model2mul, type = "LRT")


###################################################
### code chunk number 34: CNVassoc_vignette.Rnw:559-560
###################################################
coef(summary(model1mul, ref = 2))


###################################################
### code chunk number 35: CNVassoc_vignette.Rnw:569-571
###################################################
model2add <- CNVassoc(casco ~ CNV.2, data = dataMLPA, model = "add")
model2add


###################################################
### code chunk number 36: CNVassoc_vignette.Rnw:581-582
###################################################
summary(model2add)


###################################################
### code chunk number 37: CNVassoc_vignette.Rnw:595-596
###################################################
anova(model2mul, model2add)


###################################################
### code chunk number 38: CNVassoc_vignette.Rnw:643-646
###################################################
data(NeveData)
intensities <- NeveData$data
pheno <- NeveData$pheno


###################################################
### code chunk number 39: CNVassoc_vignette.Rnw:686-687
###################################################
data(NeveCalled)


###################################################
### code chunk number 40: CNVassoc_vignette.Rnw:696-697
###################################################
probs <- getProbs(NeveCalled)


###################################################
### code chunk number 41: CNVassoc_vignette.Rnw:704-705
###################################################
probs[1:5, 1:7]


###################################################
### code chunk number 42: CNVassoc_vignette.Rnw:726-727
###################################################
data(NeveRegions)


###################################################
### code chunk number 43: CNVassoc_vignette.Rnw:734-735
###################################################
probsRegions <- getProbsRegions(probs, NeveRegions, intensities)


###################################################
### code chunk number 44: CNVassoc_vignette.Rnw:743-745
###################################################
pvals <- multiCNVassoc(probsRegions, formula = "pheno ~ CNV", model = "mult", 
num.copies = 0:2, cnv.tol = 0.01)


###################################################
### code chunk number 45: CNVassoc_vignette.Rnw:755-757
###################################################
pvalsBH <- getPvalBH(pvals)
head(pvalsBH)


###################################################
### code chunk number 46: CNVassoc_vignette.Rnw:764-765
###################################################
cumsum(table(cut(pvalsBH[, 2], c(-Inf, 1e-5, 1e-4, 1e-3, 1e-2, 0.05))))


###################################################
### code chunk number 47: CNVassoc_vignette.Rnw:793-796
###################################################
data(SNPTEST)
dim(cases)
dim(controls)


###################################################
### code chunk number 48: CNVassoc_vignette.Rnw:800-801
###################################################
cases[1:10,1:11]


###################################################
### code chunk number 49: CNVassoc_vignette.Rnw:804-805
###################################################
cases[1:10,1:11]


###################################################
### code chunk number 50: CNVassoc_vignette.Rnw:809-810
###################################################
controls[1:10,1:11]


###################################################
### code chunk number 51: CNVassoc_vignette.Rnw:813-814
###################################################
controls[1:10,1:11]


###################################################
### code chunk number 52: CNVassoc_vignette.Rnw:842-850
###################################################
nSNP <- nrow(cases)
probs <- lapply(1:nSNP, function(i) {
  snpi.cases <- matrix(as.double(cases[i, 6:ncol(cases)]), ncol = 3,
byrow = TRUE)
  snpi.controls <- matrix(as.double(controls[i, 6:ncol(controls)]),
ncol = 3, byrow = TRUE)
  return(rbind(snpi.cases, snpi.controls))
})


###################################################
### code chunk number 53: CNVassoc_vignette.Rnw:859-860
###################################################
casecon <- rep(1:0, c(500, 500))


###################################################
### code chunk number 54: CNVassoc_vignette.Rnw:865-867
###################################################
pvals <- multiCNVassoc(probs, formula = "casecon~CNV", model = "add",
num.copies = 0:2, cnv.tol = 0.001)


###################################################
### code chunk number 55: CNVassoc_vignette.Rnw:871-873
###################################################
pvalsBH <- getPvalBH(pvals)
head(pvalsBH)


###################################################
### code chunk number 56: CNVassoc_vignette.Rnw:877-878
###################################################
table(cut(pvalsBH[, 2], c(-Inf, 1e-3, 1e-2, 0.05, 0.1, Inf)))


###################################################
### code chunk number 57: CNVassoc_vignette.Rnw:915-922
###################################################
set.seed(123456)
rr <- 1.7
incid0 <- 0.12
lambda <- c(incid0, incid0 * rr, incid0 * rr^2)
dsim <- simCNVdataPois(n = 4000, mu.surrog = 0:2, sd.surrog = rep(0.4, 3),
w = c(0.25, 0.5, 0.25), lambda = lambda)
head(dsim)


###################################################
### code chunk number 58: CNVassoc_vignette.Rnw:934-936
###################################################
CNV <- cnv(dsim$surrog, mix = "mclust")
CNV


###################################################
### code chunk number 59: CNVassoc_vignette.Rnw:941-943
###################################################
fit <- CNVassoc(resp ~ CNV, data = dsim, family = "poisson", model = "add")
coef(summary(fit))


###################################################
### code chunk number 60: CNVassoc_vignette.Rnw:950-955
###################################################
fit.gold <- glm(resp ~ cnv, data = dsim, family = "poisson")
table.gold <- c(exp(c(coef(fit.gold)[2], confint(fit.gold)[2,])), 
coef(summary(fit.gold))[2,4])
names(table.gold) <- c("RR", "lower", "upper", "p-value")
table.gold


###################################################
### code chunk number 61: CNVassoc_vignette.Rnw:965-970
###################################################
fit.naive <- glm(resp ~ CNV, data = dsim, family = "poisson")
table.naive <- c(exp(c(coef(fit.naive)[2], confint(fit.naive)[2,])), 
coef(summary(fit.naive))[2,4])
names(table.naive) <- c("RR", "lower", "upper", "p-value")
table.naive


###################################################
### code chunk number 62: CNVassoc_vignette.Rnw:974-980
###################################################
taula <- rbind(
table.gold[1:3],
coef(summary(fit))[,c("RR","lower.lim","upper.lim")],
table.naive[1:3])
rownames(taula)<-c("Gold","LC","Naive")
xtable(taula,"Comparison of RR estimated by the gold standard model, a latent class model (LC) and naive approach",label="table-compare")


###################################################
### code chunk number 63: CNVassoc_vignette.Rnw:1016-1030
###################################################
set.seed(123456)
n <- 5000
w <- c(0.25, 0.5, 0.25)
mu.surrog <- 0:2
sd.surrog <- rep(0.4, 3)
hr <- 1.5
incid0 <- 0.05
lambda <- c(incid0, incid0 * hr, incid0 * hr^2)
shape <- 1
scale <- lambda^(-1/shape)
perc.obs <- 0.10
time.cens <- qweibull(perc.obs, mean(shape), mean(scale))
dsim <- simCNVdataWeibull(n, mu.surrog, sd.surrog, w, lambda, shape, time.cens)
head(dsim)


###################################################
### code chunk number 64: CNVassoc_vignette.Rnw:1044-1046
###################################################
CNV<-cnv(dsim$surrog,mix="mclust")
CNV


###################################################
### code chunk number 65: CNVassoc_vignette.Rnw:1053-1055
###################################################
fit<-CNVassoc(Surv(resp,cens)~CNV, data=dsim, family="weibull", model="add")
coef(summary(fit))


