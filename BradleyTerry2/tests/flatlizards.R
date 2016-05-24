options(digits = 4) ## only applies to this file
library(BradleyTerry2)
data(flatlizards, package = "BradleyTerry2")
##
##  Fit the standard Bradley-Terry model, using the bias-reduced
##  maximum likelihood method:
##
attach(flatlizards)
result <- rep(1, nrow(contests))
BTmodel <- BTm(result, winner, loser, br = TRUE, data = contests)
summary(BTmodel)
##
##  That's fairly useless, though, because of the rather small
##  amount of data on each lizard.  And really the scientific
##  interest is not in the abilities of these particular 77
##  lizards, but in the relationship between ability and the
##  measured predictor variables.
##
##  So next fit (by maximum likelihood) a "structured" B-T model in
##  which abilities are determined by a linear predictor.
##
##  This reproduces results reported in Table 1 of Whiting et al. (2006):
##
Whiting.model <- BTm(result, winner, loser, ~ throat.PC1[..] + throat.PC3[..] +
                     head.length[..] + SVL[..], family = binomial,
                     data = list(contests, predictors))
summary(Whiting.model)
##
##  Equivalently, fit the same model using glmmPQL:
##
Whiting.model <- BTm(result, winner, loser, ~ throat.PC1[..] + throat.PC3[..] +
                     head.length[..] + SVL[..] + (1|..), sigma = 0,
                     sigma.fixed = TRUE, data = list(contests, predictors))
summary(Whiting.model)
##
##  But that analysis assumes that the linear predictor formula for
##  abilities is _perfect_, i.e., that there is no error in the linear
##  predictor.  This will always be unrealistic.
##
##  So now fit the same predictor but with a normally distributed error
##  term --- a generalized linear mixed model --- by using the BTm
##  function instead of glm.
##
Whiting.model2 <- BTm(result, winner, loser, ~ throat.PC1[..] + throat.PC3[..] +
                      head.length[..] + SVL[..] + (1|..),
                      data = list(contests, predictors), trace = TRUE)
summary(Whiting.model2)
##
##  The estimated coefficients (of throat.PC1, throat.PC3,
##  head.length and SVL are not changed substantially by
##  the recognition of an error term in the model; but the estimated
##  standard errors are larger, as expected.  The main conclusions from
##  Whiting et al. (2006) are unaffected.
##
##  With the normally distributed random error included, it is perhaps
##  at least as natural to use probit rather than logit as the link
##  function:
##
Whiting.model3 <- BTm(result, winner, loser, ~ throat.PC1[..] + throat.PC3[..] +
                      head.length[..] + SVL[..] + (1|..),
                      family = binomial(link = "probit"),
                      data = list(contests, predictors), trace = TRUE)
summary(Whiting.model3)
BTabilities(Whiting.model3)
residuals(Whiting.model3, "grouped")
##  Note the "separate" attribute here, identifying two lizards with
##  missing values of at least one predictor variable
##
##  Modulo the usual scale change between logit and probit, the results
##  are (as expected) very similar to Whiting.model2.
