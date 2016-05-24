## ---- message=FALSE, warning=FALSE---------------------------------------
library(simr)

## ----options, echo=FALSE, message=FALSE----------------------------------
simrOptions(progress=FALSE)

## ------------------------------------------------------------------------
cbpp$obs <- 1:nrow(cbpp)
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd) + (1|obs), data=cbpp,
    family=binomial)
summary(gm1)$coef

## ------------------------------------------------------------------------
doTest(gm1, fixed("period", "lr"))

## ------------------------------------------------------------------------
doTest(gm1, fixed("period2", "z"))

## ------------------------------------------------------------------------
gm2 <- glmer(cbind(incidence, size - incidence) ~ period + size + (1 | herd), data=cbpp,
    family=binomial)
doTest(gm2, fixed("size", "z"))

## ------------------------------------------------------------------------
fixef(gm2)["size"] <- 0.05
powerSim(gm2, fixed("size", "z"), nsim=50)

## ------------------------------------------------------------------------
fm1 <- lmer(angle ~ recipe * temp + (1|recipe:replicate), data=cake, REML=FALSE)

## ------------------------------------------------------------------------
doTest(fm1, fcompare(~ recipe + temp))

## ------------------------------------------------------------------------
fm2 <- lmer(angle ~ recipe + poly(temp, 2) + (1|recipe:replicate), data=cake, REML=FALSE)
summary(fm2)$coef
doTest(fm2, fcompare(~ recipe + temp))

## ------------------------------------------------------------------------
data(budworm, package="pbkrtest")
bw1 <- glm(cbind(ndead, ntotal-ndead) ~ dose*sex, family="binomial", data=budworm)
summary(bw1)$coef

## ------------------------------------------------------------------------
doTest(bw1, compare(. ~ dose + sex))

## ------------------------------------------------------------------------
doTest(bw1, fixed("dose:sexmale", "z"))

## ------------------------------------------------------------------------
re1 <- lmer(Yield ~ 1|Batch, data=Dyestuff)
doTest(re1, random())

## ------------------------------------------------------------------------
fm1 <- lmer(Reaction ~ Days + (Days | Subject), data=sleepstudy)

## ------------------------------------------------------------------------
doTest(fm1, compare( ~ Days + (1 | Subject)))

## ---- eval=FALSE---------------------------------------------------------
#  doTest(fm1, compare( ~ Days + (1 | Subject), "pb"))

## ---- eval=FALSE---------------------------------------------------------
#  doTest(fm1, rcompare( ~ (1 | Subject), "pb"))

## ------------------------------------------------------------------------
binFit <- glm(formula = cbind(z, 10 - z) ~ x + g, family = binomial, data = simdata)

poisSim <- glm(formula = z ~ x + g, family = poisson, data = simdata)
coef(poisSim)[1:2] <- c(1, -0.05)

powerSim(binFit, sim=poisSim, nsim=50, seed=1)

## ------------------------------------------------------------------------
ps <- lastResult()
ps$errors

