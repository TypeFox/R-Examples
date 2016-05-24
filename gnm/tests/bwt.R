library(MASS)
example(birthwt)
library(nnet)
bwt.mu <- multinom(low ~ ., data = bwt)

## Equivalent using gnm - include unestimable main effects in model so
## that interactions with low0 automatically set to zero, else could use
## 'constrain' argument.
library(gnm)
bwtLong <- expandCategorical(bwt, "low", group = FALSE)
bwt.po <- gnm(count ~  low*(. - id), eliminate = id, data = bwtLong, family =
              "poisson")

coef(bwt.po)

summary(bwt.po)

anova(bwt.po)

drop1(bwt.po)

bwt.po <- gnm(count ~  low*age - id, eliminate = id, data = bwtLong, family =
              "poisson")

add1(bwt.po, formula(terms(count~low*( . -id), data = bwtLong)))

bwt.po <- gnm(count ~  . - id, eliminate = id, data = bwtLong, family =
              "poisson")


