### R code from vignette source 'phia.Rnw'

###################################################
### code chunk number 1: phia.Rnw:13-20
###################################################
options(useFancyQuotes = FALSE, width=85)
car.citation <- citation("car")
car.citation$key <- "Fox-CAR2011"
nlme.citation <- citation("nlme")[1]
nlme.citation$key <- "Pinheiro-nlme"
lme4.citation <- citation("lme4")[1]
lme4.citation$key <- "Bates-lme4"


###################################################
### code chunk number 2: phia.Rnw:24-27
###################################################
toBibtex(car.citation)
toBibtex(nlme.citation)
toBibtex(lme4.citation)


###################################################
### code chunk number 3: phia.Rnw:238-240
###################################################
library(phia)
some(Boik)


###################################################
### code chunk number 4: modboikresid
###################################################
mod.boik <- lm(edr ~ therapy*medication, data=Boik)
par(mfcol=c(1,2))
plot(mod.boik, 1:2) # Plot diagnostics for the model
Anova(mod.boik)


###################################################
### code chunk number 5: phia.Rnw:259-260
###################################################
mod.boik <- lm(edr ~ therapy*medication, data=Boik)
par(mfcol=c(1,2))
plot(mod.boik, 1:2) # Plot diagnostics for the model
Anova(mod.boik)


###################################################
### code chunk number 6: phia.Rnw:295-296
###################################################
(boik.means <- interactionMeans(mod.boik))


###################################################
### code chunk number 7: phia.Rnw:309-310
###################################################
interactionMeans(mod.boik, factors="therapy")


###################################################
### code chunk number 8: phia.Rnw:329-330
###################################################
plot(boik.means)


###################################################
### code chunk number 9: phia.Rnw:345-346
###################################################
plot(boik.means, multiple=FALSE) # Not printed in this paper


###################################################
### code chunk number 10: phia.Rnw:416-417
###################################################
testInteractions(mod.boik, fixed="therapy", across="medication")


###################################################
### code chunk number 11: phia.Rnw:473-476
###################################################
boik.mtable <- xtabs(boik.means$"adjusted mean" ~ therapy+medication, boik.means)
boik.mtable <- addmargins(boik.mtable, FUN=mean, quiet=TRUE)
print(boik.mtable, digits=4)


###################################################
### code chunk number 12: phia.Rnw:484-488
###################################################
boik.resid <- boik.mtable - boik.mtable[4,5] # Subtract the mean
boik.resid <- sweep(boik.resid, 1, boik.resid[,5]) # Subtract row means
boik.resid <- sweep(boik.resid, 2, boik.resid[4,]) # Subtract column means
print(boik.resid, digits=4)


###################################################
### code chunk number 13: phia.Rnw:495-496
###################################################
testInteractions(mod.boik,residual=c("therapy","medication"))


###################################################
### code chunk number 14: residualsplot
###################################################
matplot(t(boik.resid[-4,-5]), type="b", xaxt="n", ylab="Interaction residuals")
axis(1, at=1:4, labels=levels(Boik$medication))


###################################################
### code chunk number 15: phia.Rnw:512-513
###################################################
matplot(t(boik.resid[-4,-5]), type="b", xaxt="n", ylab="Interaction residuals")
axis(1, at=1:4, labels=levels(Boik$medication))


###################################################
### code chunk number 16: phia.Rnw:580-581
###################################################
testInteractions(mod.boik, pairwise="therapy", across="medication")


###################################################
### code chunk number 17: phia.Rnw:598-599
###################################################
testInteractions(mod.boik)


###################################################
### code chunk number 18: phia.Rnw:635-642
###################################################
(custom.contr <- contrastCoefficients(
therapy ~ control - (T1 + T2)/2,          # Control vs. both therapies
therapy ~ T1 - T2,                        # Therapy T1 vs. T2
medication ~ placebo - (D1 + D2 + D3)/3,  # Placebo vs. all doses
medication ~ D1 - D3,                     # Min. dose vs. max
medication ~ D2 - (D1 + D2 + D3)/3,       # Med. dose vs. average
data=Boik, normalize=TRUE))               # Normalize to homogeinize the scale


###################################################
### code chunk number 19: phia.Rnw:649-652
###################################################
names(custom.contr$therapy) <- c("cntrl.vs.all", "T1.vs.T2")
names(custom.contr$medication) <- c("plcb.vs.all", "D1.vs.D3", "D2.vs.avg")
testInteractions(mod.boik,custom=custom.contr)


###################################################
### code chunk number 20: phia.Rnw:690-691
###################################################
some(OBrienKaiser, 6)


###################################################
### code chunk number 21: phia.Rnw:703-704
###################################################
(idata <- expand.grid(hour=ordered(1:5), phase=c("pre", "post","fup")))


###################################################
### code chunk number 22: phia.Rnw:710-711
###################################################
addmargins(table(OBrienKaiser[c("gender","treatment")]))


###################################################
### code chunk number 23: phia.Rnw:717-721
###################################################
mod.ok <- lm(cbind(pre.1, pre.2, pre.3, pre.4, pre.5,
    post.1, post.2, post.3, post.4, post.5,
    fup.1, fup.2, fup.3, fup.4, fup.5) ~ treatment*gender,
    data=OBrienKaiser)


###################################################
### code chunk number 24: phia.Rnw:730-731
###################################################
Anova(mod.ok, idata=idata, idesign=~phase*hour, type=3)


###################################################
### code chunk number 25: okplot
###################################################
ok.means <- interactionMeans(mod.ok, c("treatment","phase"), idata=idata)
plot(ok.means, errorbar="ci95")


###################################################
### code chunk number 26: phia.Rnw:766-767
###################################################
ok.means <- interactionMeans(mod.ok, c("treatment","phase"), idata=idata)
plot(ok.means, errorbar="ci95")


###################################################
### code chunk number 27: phia.Rnw:785-786
###################################################
testInteractions(mod.ok, pairwise=c("treatment", "phase"), idata=idata)


###################################################
### code chunk number 28: phia.Rnw:861-862
###################################################
str(Prestige)


###################################################
### code chunk number 29: phia.Rnw:875-877
###################################################
mod.prestige <- lm(prestige ~ (log2(income)+education)*type, data=Prestige)
Anova(mod.prestige)


###################################################
### code chunk number 30: prestigeplot
###################################################
plot(interactionMeans(mod.prestige, atx="type", slope="log2(income)"))
testInteractions(mod.prestige, pairwise="type", slope="log2(income)")


###################################################
### code chunk number 31: phia.Rnw:896-897
###################################################
plot(interactionMeans(mod.prestige, atx="type", slope="log2(income)"))
testInteractions(mod.prestige, pairwise="type", slope="log2(income)")


###################################################
### code chunk number 32: phia.Rnw:915-916
###################################################
table(Prestige$type) # Frequencies of occupation types


###################################################
### code chunk number 33: phia.Rnw:940-943
###################################################
mod.prestige2 <- update(mod.prestige, formula=.~.+(log2(income):education)*type,
subset=(Prestige$type!="wc"))
Anova(mod.prestige2)


###################################################
### code chunk number 34: phia.Rnw:951-952
###################################################
testInteractions(mod.prestige2, pairwise="type", slope="log2(income)")


###################################################
### code chunk number 35: phia.Rnw:966-970
###################################################
# Look quantiles of the model frame (a subset of the original data)
quantile(model.frame(mod.prestige2)$education)
testInteractions(mod.prestige2, pairwise="type", slope="log2(income)",
covariates=c(education=14))


###################################################
### code chunk number 36: phia.Rnw:1031-1032
###################################################
str(AMSsurvey)


###################################################
### code chunk number 37: phia.Rnw:1055-1056
###################################################
ftable(xtabs(count ~ sex + citizen + type, AMSsurvey))


###################################################
### code chunk number 38: phia.Rnw:1066-1067
###################################################
mod.ams <- glm(count ~ type*(sex+citizen), family=poisson, data=AMSsurvey)


###################################################
### code chunk number 39: phia.Rnw:1073-1074
###################################################
Anova(mod.ams)


###################################################
### code chunk number 40: amsplot
###################################################
ams.means <- interactionMeans(mod.ams)
plot(ams.means, atx="type", traces=c("sex","citizen"))


###################################################
### code chunk number 41: phia.Rnw:1092-1093
###################################################
ams.means <- interactionMeans(mod.ams)
plot(ams.means, atx="type", traces=c("sex","citizen"))


###################################################
### code chunk number 42: phia.Rnw:1125-1127
###################################################
testInteractions(mod.ams, pairwise=c("type","sex")) # test type:sex
testInteractions(mod.ams, pairwise=c("type","citizen")) #test type:citizen


###################################################
### code chunk number 43: phia.Rnw:1176-1177 (eval = FALSE)
###################################################
## dm.de <- family(model)$mu.eta


###################################################
### code chunk number 44: phia.Rnw:1302-1308
###################################################
Snijders <- nlme::bdf[c("langPRET","langPOST",   # Outcomes
    "pupilNR", "IQ.ver.cen", "ses", "sex",       # Student-related variables
    "schoolNR", "schoolSES", "avg.IQ.ver.cen")]  # School-related variables
Snijders$sex <- factor(Snijders$sex, labels=c("F","M"))
names(Snijders) <-
    c("score.1","score.2","student","IQ","SES","sex","school","avgSES","avgIQ")


###################################################
### code chunk number 45: phia.Rnw:1315-1317
###################################################
Snijders$IQ2.pos <- with(Snijders, (IQ > 0)*IQ^2)
Snijders$IQ2.neg <- with(Snijders, (IQ < 0)*IQ^2)


###################################################
### code chunk number 46: phia.Rnw:1328-1329
###################################################
library(lme4)


###################################################
### code chunk number 47: phia.Rnw:1331-1335 (eval = FALSE)
###################################################
## form1 <- score.2 ~ 
##     IQ * SES + IQ2.pos + IQ2.neg + sex + avgIQ * avgSES +    # Fixed part
##     (IQ | school)                                            # Random part
## mod.snijders.1 <- lmer(form1, data=Snijders)


###################################################
### code chunk number 48: phia.Rnw:1347-1349 (eval = FALSE)
###################################################
## form2.1 <- update(form1, (score.1+score.2)/2~.)
## form2.2 <- update(form1, (score.2-score.1)~.)


###################################################
### code chunk number 49: phia.Rnw:1368-1374
###################################################
Snijders.long <- reshape(Snijders, direction="long", idvar="student",
    varying=list(c("score.1","score.2")), v.names="score", timevar="repetition")
# The within-subjects factor must be coded as a factor
Snijders.long$repetition <- as.factor(Snijders.long$repetition)
# See the variables of the long data frame
str(Snijders.long)


###################################################
### code chunk number 50: phia.Rnw:1390-1397
###################################################
form3 <- score ~ 
    repetition * (IQ * SES + IQ2.pos + IQ2.neg + sex) +   # Student-related
    avgIQ * avgSES +                                      # School-related
    (IQ | school) + (1 | student)                         # Random part
mod.snijders.3 <- lmer(form3, data=Snijders.long)
# See the main parameters of the model (ommit correlations table)
print(mod.snijders.3, correlation=FALSE)


###################################################
### code chunk number 51: phia.Rnw:1428-1429
###################################################
Anova(mod.snijders.3)


###################################################
### code chunk number 52: phia.Rnw:1442-1448
###################################################
# Cell means
interactionMeans(mod.snijders.3)
# Simple effect of sex at each repetition
testInteractions(mod.snijders.3, fixed="repetition", across="sex")
# Pairwise interactions (default)
testInteractions(mod.snijders.3)


###################################################
### code chunk number 53: phia.Rnw:1481-1482
###################################################
graphics.off()


