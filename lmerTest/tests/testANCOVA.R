require(lmerTest)

testType1 <- TRUE

if(testType1){
load(system.file("testdata", "tree.RData", package="lmerTest"))


lmer1 <- lmer(increase ~ height*treat + (1|block), data = tree)
step1 <- step(lmer1)

an1 <- anova(lmer1, type=1)
an1.lme4 <- anova(lmer1, ddf="lme4")

stopifnot(all.equal(an1[, "Sum Sq"], an1.lme4[,"Sum Sq"]),
          all.equal(an1[, "Mean Sq"], an1.lme4[,"Mean Sq"]),
          all.equal(an1[, "F.value"], an1.lme4[,"F value"]))


lmer2 <- lmer(increase ~ treat-1+ height:treat - height + (1|block), data = tree)

an2.3 <- anova(lmer2)
an2.1 <- anova(lmer2, type=1)

TOL <- 1e-5
## check type3
stopifnot(all.equal(an2.3[, "F.value"], c(28.6195, 40.7608), tol= TOL),
          all.equal(an2.3[, "DenDF"], c(27, 27), tol= TOL))
##check type1
stopifnot(all.equal(an2.1[, "F.value"], c(787.736, 40.7608), tol= TOL),
          all.equal(an2.1[, "DenDF"], c(27, 27), tol= TOL))

## check for the model with one factor
lmer4 <- lmer(increase ~ treat - 1+ (1|block), data = tree)

an4.1 <- anova(lmer4, type=1)
stopifnot(all.equal(an4.1[, "F.value"], c(80.4380), tol= TOL))



lmer5 <- lmer(increase ~ treat - 1 + height:treat+ (1|block), data = tree)

an5.1 <- anova(lmer5, type=1)

stopifnot(all.equal(an2.1, an5.1))

## check for one covariate
lmer6 <- lmer(increase ~  height - 1+ (1|block), data = tree)

an6.1 <- anova(lmer6, type=1)

TOL <- 1e-3
stopifnot(all.equal(an6.1[,"F.value"], 347.374, tol=TOL),
          all.equal(an6.1[,"DenDF"], 10.22, tol=TOL))

## check for multiple factors
m <- lmer(Coloursaturation ~ TVset*Picture - 1 +
            (1|Assessor), data=TVbo)
an <- anova(m)
an.1 <- anova(m, type = 1)

TOL <- 1e-6
TOL2 <- 1e-3
stopifnot(all.equal(an[, "F.value"], c(84.5327, 2.85346, 1.75537), tol = TOL),
          all.equal(an.1[, "F.value"], c(615.423, 2.85346, 1.75537), tol = TOL),
          all.equal(an.1[, "DenDF"], c(16.17, 173, 173), tol = TOL2))

m.carrots <- lmer(Preference ~ sens2 + Homesize - 1
                  +(1+sens2|Consumer), data=carrots)

an.1 <- anova(m.carrots, type=1)
an.3 <- anova(m.carrots)
an.lme4 <- anova(m.carrots, ddf = "lme4")

TOL <- 1e-5
stopifnot(all.equal(an.1[, "F.value"], c(56.5394, 4169.87), tol = TOL),
          all.equal(an.3[, "F.value"], c(54.8206, 4169.87), tol = TOL))


# the one from lme4 is different
tools::assertError(stopifnot(all.equal(an.lme4[, "F value"], an.1, tol = TOL)))


## checking lsmeans
lmer4.noint <- lmer(increase ~ treat - 1 + (1|block), data = tree)
lmer4 <- lmer(increase ~ treat + (1|block), data = tree)
stopifnot(all.equal(lsmeans(lmer4), lsmeans(lmer4.noint)))

lsm <- lsmeans(lmer4)

TOL <- 1e-5
stopifnot(all.equal(lsm$lsmeans.table[, "Estimate"], 
                    as.numeric(tapply(tree$increase, tree$treat, mean)), 
                    check.names = FALSE, check.attributes = FALSE,tol = TOL))

}