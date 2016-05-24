library("Countr")
         # object <- readRDS("McShane_Wei_results_boot.RDS")
fn <- system.file("extdata", "McShane_Wei_results_boot.RDS", package = "Countr")
object <- readRDS(fn)
data <- object$data


##-------------------- test coef renewal
print("testing coef renewal ....")
print(coef(object))
Sys.sleep(10)

##-------------------- test vcov renewal
print("testing vcov renewal ....")
print(vcov(object))
Sys.sleep(10)

## --------------------- test residuals renewal
print("testing residuals ...")
print("**** Pearson residuals: rescaled by sd ****")
print(head(residuals(object, "pearson")))
print("**** response residuals: non rescaled ****")
print(head(residuals(object, "response")))
Sys.sleep(10)

## ------------------- test fitted value renewal
print("testing fitted values (mean response) ...")
print(head(fitted(object)))
Sys.sleep(10)

##------------------- test coef std errors
print("testing coefficients standard errors ...")
asym <- se.coef(object, , "asymptotic")
boot <- se.coef(object, , "boot")
print(cbind(asym, boot))
Sys.sleep(10)

##------------------ test confint
print("testing coefficients ci ...")
asym <- confint(object, type = "asymptotic")
boot <- confint(object, type = "boot", bootType = "norm")
print(list(asym = asym, boot = boot))
Sys.sleep(10)

## -------------------- summary renewal
print(summary(object))
Sys.sleep(10)

##----------------------- print method
print(object)
Sys.sleep(10)

##----------------------- loglik, nobs, AIC, BIC method
print("testing loglik, nobs, AIC, BIC ....")
print(c(loglik = as.numeric(logLik(object)), nobs = nobs(object),
        AIC = AIC(object), BIC = BIC(object)))
Sys.sleep(10)

##-------------------------- predict method
print("testing prediction method ...")
## old data
predOld.response <- predict(object, type = "response", se.fit = TRUE)
predOld.prob <- predict(object, type = "prob", se.fit = TRUE)

## newData (extracted from old Data)
newData <- head(data)
predNew.response <- predict(object, newdata = newData,
                            type = "response", se.fit = TRUE)
predNew.prob <- predict(object, newdata = newData,
                        type = "prob", se.fit = TRUE)

print(cbind(head(predOld.response$values),
            head(predOld.response$se$scale),
            head(predOld.response$se$shape),
            predNew.response$values,
            predNew.response$se$scale,
            predNew.response$se$shape)
      )
Sys.sleep(10)

print(cbind(head(predOld.prob$values),
            head(predOld.prob$se$scale),
            head(predOld.prob$se$shape),
            predNew.prob$values,
            predNew.prob$se$scale,
            predNew.prob$se$shape)
      )
