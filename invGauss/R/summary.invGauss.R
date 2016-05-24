"summary.invGauss" <-
function(object, covariance = FALSE, ...)
{
##
## EXTRACTS AND COMPUTES MOST IMPORTANT ELEMENTS OF THE invGauss FIT
##
.call <- object$call
.coef <- object$coefficients
.var <- object$cov.unscaled
.loglik <- object$loglik
.AIC <- object$AIC

.SE <-  sqrt(diag(.var))
.Z <- .coef/.SE
.pval <- 2*pnorm(abs(.Z), lower.tail = F) # two-sided test

.tab <- cbind(Estimate = .coef, "Std. Error" = .SE, "Z-score" = .Z, "Pr(>|z|)" = .pval)
.tab["tau", "Pr(>|z|)"] <- NA # Probably not a good idea always to trust a Wald test


.ut <- list(call = .call, coefficients = .tab, cov.unscaled = .var, loglik = .loglik, AIC = .AIC)
class(.ut) <- "summary.invGauss"
return(.ut)
}

