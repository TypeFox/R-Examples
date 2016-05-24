date()
data(armd, package = "nlmeU")
library(nlme)
lm3.form <- formula(visual ~ visual0 + time + treat.f)
fm16.5ml <- lme(lm3.form, random = list(subject = pdDiag(~time)) ,
               weights = varPower(form = ~ time),  data = armd,
               method="ML")


df1 <- subset(armd, subject %in% "1") # Data for subject "1"
detach(package:nlme)

library(nlmeU)
logLik1(fm16.5ml, df1)


packageVersion("nlme")
sessionInfo()
detach(package:nlmeU)
