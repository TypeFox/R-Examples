date()
data(armd, package = "nlmeU")
library(nlme)
lm3.form <- formula(visual ~ visual0 + time + treat.f)
fm16.5ml <- lme(lm3.form, random = list(subject = pdDiag(~time)) ,
               weights = varPower(form = ~ time),  data = armd,
               method="ML")
detach(package:nlme)

library(nlmeU)
simY <- simulateY(fm16.5ml, nsim = 10, seed = 5342)
str(simY)

packageVersion("nlme")
sessionInfo()
detach(package:nlmeU)