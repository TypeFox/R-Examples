

#################################################################################

library(msme)

# library(msme, lib.loc="lib")

library(MASS)

data(medpar)

denom <- rep(1:5, each=299, times=1)*100   # m : binomial denominator w medpar
oset <- rep(1:5, each=299, times=1)*100    # offset Poisson, NB, offset w medpar
loset <- log(oset)                         # log of oset

## POISSON -------------------------------------------------

glm.poi <- glm(los ~ hmo + white,
               family = "poisson",
               data = medpar)

ml.poi <- ml_glm(los ~ hmo + white,
                 family = "poisson",
                 link = "log",
                 data = medpar)

ml.poi  
glm.poi

summary(ml.poi)
summary(glm.poi)

## RATE POISSON

glm.rpoi <- glm(los ~ hmo + white + offset(loset),
                family = poisson,
                data = medpar)
ml.rpoi <- ml_glm(los ~ hmo + white,
                  family = "poisson", link = "log",
                  offset = loset,
                  data = medpar)

ml.rpoi
glm.rpoi

summary(ml.rpoi)
summary(glm.rpoi)

## IDENTITY POISSON ----------------------------old-----------

glm.ipoi <- glm(los ~ hmo + white,
                family = poisson(link=identity),
                data = medpar)

ml.ipoi <- ml_glm(los ~ hmo + white,
                  family = "poisson", link = "identity",
                  data = medpar)

glm.ipoi
ml.ipoi

summary(glm.ipoi)
summary(ml.ipoi)

##  BERNOULLI LOGIT --------------------------------------------

glm.logit <- glm(died ~ hmo + white,
                 family = "binomial",
                 data = medpar)
ml.logit <- ml_glm(died ~ hmo + white,
                  data = medpar,
                  family = "bernoulli",
                  link = "logit1")

ml.logit
glm.logit

summary(ml.logit)
summary(glm.logit)

## BERNOULLI PROBIT --------------------------------------------

ml.probit <- ml_glm(died ~ hmo + white,
                    family = "bernoulli",
                    link = "probit1",
                    data = medpar)
glm.probit <- glm(died ~ hmo + white,
                  family = binomial(link=probit),
                  data = medpar)


ml.probit
glm.probit

summary(ml.probit)
summary(glm.probit)

## BERNOULLI CLOGLOG --------------------------------------------

ml.cloglog <- ml_glm(died ~ hmo + white,
                     family = "bernoulli",
                     link = "cloglog1",
                     data = medpar,
                     verbose = 1)

glm.cloglog <- glm(died ~ hmo + white,
                   family = binomial(link=cloglog),
                   data = medpar)

ml.cloglog
glm.cloglog

summary(ml.cloglog)
summary(glm.cloglog)

