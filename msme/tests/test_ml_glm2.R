library(msme)

# library(msme, lib.loc="lib")

library(MASS)

data(medpar)
data(ufc)

denom <- rep(1:5, each=299, times=1)*100   # m : binomial denominator w medpar
oset <- rep(1:5, each=299, times=1)*100    # offset Poisson, NB, offset w medpar
loset <- log(oset)                         # log of oset

## Normal

ufc <- na.omit(ufc)
ml.g <- ml_glm2(height.m ~ dbh.cm,
                formula2 = ~1,
                data = ufc,
                family = "normal",
                mean.link = "identity",
                scale.link = "log_s")

lm.g <- lm(height.m ~ dbh.cm,
                data = ufc)
                
ml.g
lm.g

summary(ml.g)
summary(lm.g)

## NEGATIVE BINOMIAL (NB2) ------------------------------------------

glm.nb2 <- glm.nb(los ~ hmo + white,
                  data = medpar)

ml.nb2 <- ml_glm2(los ~ hmo + white,
                    formula2 = ~1,
                    data = medpar,
                    family = "negBinomial",
                    mean.link = "log",
                    scale.link = "log_s")

glm.nb2
ml.nb2

summary(glm.nb2)
summary(ml.nb2)

## RATE NEGATIVE BINOMIAL (NB2)--------------------------------

glm.rnb2 <- glm.nb(los ~ hmo + white + offset(loset),
                   data = medpar)

ml.rnb2 <- ml_glm2(los ~ hmo + white,
                    formula2 = ~1,
                    data = medpar,
                  offset=loset,
                    family = "negBinomial",
                    mean.link = "log",
                    scale.link = "inverse_s")

ml.rnb2
glm.rnb2

summary(glm.rnb2)
summary(ml.rnb2)



## GAMMA CANONICAL--------------------------------------
#irls.gam <- irls(los ~ hmo + white,
#                 family = "gamma",
#                 link = "inverse",
#                 data = medpar)

#glm.gam <- glm(los ~ hmo + white,
#               family = Gamma,
#               data = medpar)

#ml.gam <- ml_glm2(los ~ hmo + white,
#                    formula2 = ~1, 
#                    data = medpar,
#                    family = "gamma",
#                    mean.link = "inverse",
#                    scale.link = "inverse_s")

#irls.gam
#glm.gam
#ml.gam

#summary(irls.gam)
#summary(glm.gam)
#summary(ml.gam)


## INVERSE GAUSSIAN CANONICAL--------------------------------------
#irls.ivg <- irls(los ~ hmo + white,
#                 family = "inv_gauss",
#                 link = "inverse2",
#                 data = medpar)
#glm.ivg <- glm(los ~ hmo + white,
#               family = inverse.gaussian,
#               data = medpar)

#irls.ivg
#glm.ivg

#summary(irls.ivg)
#summary(glm.ivg)

