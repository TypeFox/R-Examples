require(lmerTest)
testType1 <- TRUE


#if(testType1){
## from merModLmerTest
m <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject),
          data = sleepstudy)

anova(m, type=1)

if(require(pbkrtest))
  anova(m, ddf="Kenward-Roger")

## from lmerTest
m <- lmer(Informed.liking ~ Gender+Information+Product +(1|Consumer), data=ham)

if(require(pbkrtest))
  anova(m, ddf="Kenward-Roger")

## from anova methods
m.ham <- lmer(Informed.liking ~ Product*Information*Gender 
              + (1|Consumer), data = ham)

an.lmerTest <- anova(m.ham, type = 1)

an.lme4 <- anova(m.ham, ddf = "lme4")

## check that F values and so SS and MS agree in lmerTest and lme4
stopifnot(all.equal(an.lmerTest[,"F.value"], an.lme4[, "F value"]),
          all.equal(an.lmerTest[,"Sum Sq"], an.lme4[, "Sum Sq"]),
          all.equal(an.lmerTest[,"Mean Sq"], an.lme4[, "Mean Sq"]))

fm2 <- lmer(Preference ~ sens2 + I(sens1^2)  +
              (1+sens2|Consumer), data=carrots)

## from lmer
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0 + Days|Subject), sleepstudy)


an.sleep <- anova(fm1, type = 1)
an.sleep2 <- anova(fm2, type = 1)


## check with SAS
TOL <- 1e-6
stopifnot(all.equal(an.sleep[, "F.value"], 45.8530, tol = TOL))
stopifnot(all.equal(an.sleep[, "DenDF"], 17, tol = TOL))

## for model withut correlations
stopifnot(all.equal(an.sleep[, "F.value"], 45.8530, tol = TOL))
stopifnot(all.equal(an.sleep[, "DenDF"], 17, tol = TOL))

## check KR works
anova(fm1, ddf="Kenward-Roger")

# anova table the same as of class merMod
anova(fm1, ddf="lme4")

anova(fm1, fm2)

#}
