date()
library(testthat)
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
res <- Pwr(fm16.5,  L = c("treat.fActive" = 1))   # The L argument
library(testthat)
expect_that(unlist(res["F-value"]), equals(c("F-value" = 5.534487163)))
packageVersion("nlme")
sessionInfo()
detach(package:nlmeU)
