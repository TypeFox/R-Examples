### R code from vignette source 'Usinglmer.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(width=75, contrasts=c(unordered="contr.SAS",ordered="contr.poly"))
library(SASmixed)
library(lme4)
library(lattice)


###################################################
### code chunk number 2: applyClass
###################################################
sapply(Animal, data.class)
str(Animal)


###################################################
### code chunk number 3: semiconductorGrp
###################################################
Semiconductor <- within(Semiconductor, Grp <- factor(ET:Wafer))


###################################################
### code chunk number 4: contrasts (eval = FALSE)
###################################################
## options(contrasts = c(factor = "contr.SAS", ordered = "contr.poly"))


###################################################
### code chunk number 5: adg1
###################################################
print(xyplot(adg ~ Treatment | Block, AvgDailyGain, type = c("g", "p", "r"),
       xlab = "Treatment (amount of feed additive)",
       ylab = "Average daily weight gain (lb.)", aspect = "xy",
       index.cond = function(x, y) coef(lm(y ~ x))[1]))


###################################################
### code chunk number 6: adg
###################################################
## compare with output 5.1, p. 178
(fm1Adg <- lmer(adg ~ (Treatment - 1)*InitWt + (1 | Block), AvgDailyGain))
anova(fm1Adg)   # checking significance of terms
## common slope model
(fm2Adg <- lmer(adg ~ InitWt + Treatment + (1 | Block), AvgDailyGain))
anova(fm2Adg)
(fm3Adg <- lmer(adg ~ InitWt + Treatment - 1 + (1 | Block), AvgDailyGain))


###################################################
### code chunk number 7: bib1
###################################################
print(xyplot(y ~ x | Block, BIB, groups = Treatment, type = c("g", "p"),
             aspect = "xy", auto.key = list(points = TRUE, space = "right",
             lines = FALSE)))


###################################################
### code chunk number 8: bib
###################################################
## compare with Output 5.7, p. 188
(fm1BIB <- lmer(y ~ Treatment * x + (1|Block), BIB))
anova(fm1BIB)     # strong evidence of different slopes
## compare with Output 5.9, p. 193
(fm2BIB <- lmer(y ~ Treatment + x:Grp + (1|Block), BIB))
anova(fm2BIB)


###################################################
### code chunk number 9: bond
###################################################
## compare with output 1.1 on p. 6
(fm1Bond <- lmer(pressure ~ Metal + (1|Ingot), Bond))
anova(fm1Bond)


###################################################
### code chunk number 10: Cultivation
###################################################
str(Cultivation)
xtabs(~Block+Cult, Cultivation)
(fm1Cult <- lmer(drywt ~ Inoc * Cult + (1|Block) + (1|Cult), Cultivation))
anova(fm1Cult)
(fm2Cult <- lmer(drywt ~ Inoc + Cult + (1|Block) + (1|Cult), Cultivation))
anova(fm2Cult)
(fm3Cult <- lmer(drywt ~ Inoc + (1|Block) + (1|Cult), Cultivation))
anova(fm3Cult)


###################################################
### code chunk number 11: Demand
###################################################
## compare to output 3.13, p. 132
(fm1Demand <-
 lmer(log(d) ~ log(y) + log(rd) + log(rt) + log(rs) + (1|State) + (1|Year),
      Demand))


###################################################
### code chunk number 12: HR
###################################################
## linear trend in time
(fm1HR <- lmer(HR ~ Time * Drug + baseHR + (Time|Patient), HR))
anova(fm1HR)
## remove interaction
(fm3HR <- lmer(HR ~ Time + Drug + baseHR + (Time|Patient), HR))
anova(fm3HR)
## remove Drug term
(fm4HR <- lmer(HR ~ Time + baseHR + (Time|Patient), HR))
anova(fm4HR)


###################################################
### code chunk number 13: Mississippi
###################################################
## compare with output 4.1, p. 142
(fm1Miss <- lmer(y ~ 1 + (1 | influent), Mississippi))
## compare with output 4.2, p. 143
(fm1MLMiss <- lmer(y ~ 1 + (1 | influent), Mississippi, REML=FALSE))
ranef(fm1MLMiss)          # BLUP's of random effects on p. 144
ranef(fm1Miss)            # BLUP's of random effects on p. 142
VarCorr(fm1Miss)          # compare to output 4.7, p. 148
## compare to output 4.8 and 4.9, pp. 150-152
(fm2Miss <- lmer(y ~ Type + (1 | influent), Mississippi))
anova(fm2Miss)


###################################################
### code chunk number 14: Multilocation
###################################################
str(Multilocation)
### Create a Block %in% Location factor
Multilocation$Grp <- with(Multilocation, Block:Location)
(fm1Mult <- lmer(Adj ~ Location * Trt + (1|Grp), Multilocation))
anova(fm1Mult)
(fm2Mult <- lmer(Adj ~ Location + Trt + (1|Grp), Multilocation))
(fm3Mult <- lmer(Adj ~ Location       + (1|Grp), Multilocation))
(fm4Mult <- lmer(Adj ~            Trt + (1|Grp), Multilocation))
(fm5Mult <- lmer(Adj ~ 1              + (1|Grp), Multilocation))
anova(fm2Mult)
fm2MultR <- lmer(Adj ~ Trt + (Trt - 1|Location) + (1|Block), Multilocation,
                 verbose = TRUE)
## non convergence in 10000 evaluations
fm2MultR


###################################################
### code chunk number 15: PBIB
###################################################
str(PBIB)
## compare with output 1.7  pp. 24-25
(fm1PBIB <- lmer(response ~ Treatment + (1 | Block), PBIB))


###################################################
### code chunk number 16: SIMS
###################################################
str(SIMS)
## compare to output 7.4, p. 262
(fm1SIMS <- lmer(Gain ~ Pretot + (Pretot | Class), SIMS))


