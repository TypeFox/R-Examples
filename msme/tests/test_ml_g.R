

library(msme)

data(ufc)
ufc <- na.omit(ufc)

ufc.g.reg <- ml_g(height.m ~ dbh.cm, data = ufc)

ufc.g.lm <- lm(height.m ~ dbh.cm, data = ufc)

ufc.g.reg
ufc.g.lm

summary(ufc.g.reg)
summary(ufc.g.lm)

coef(ufc.g.reg)

logLik(ufc.g.reg)
logLik(ufc.g.lm)

AIC(ufc.g.reg)
AIC(ufc.g.lm)




