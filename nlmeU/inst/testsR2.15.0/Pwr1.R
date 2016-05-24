date()
## See R20.14b
data(armd, package = "nlmeU")

library(nlme)
lm3.form <- formula(visual ~ visual0 + time + treat.f) 
fm16.5 <- 
   lme(lm3.form,             
       random = list(subject = pdDiag(~time)),       
       weights = varPower(form = ~time),
       data = armd)     
formula(fm16.5)                            # Recall formula
fixef(fm16.5)
detach(package:nlme)

library(nlmeU)
Pwr(fm16.5)                                # Default call 
Pwr(fm16.5,  L = c("treat.fActive" = 1))   # The L argument
packageVersion("nlme")
sessionInfo()
detach(package:nlmeU)
