# AIC and BIC statistics following glm, glm.nb, nbinomial
# Joseph Hilbe, Modeling Count Data, Cambridge Univ Press
# version 2  13 Aug, 2014. Amend xvars line
modelfit  <- function(x)  {
obs    <- x$df.null + 1
aic    <- x$aic
xvars  <- x$df.null - x$df.residual + 1
rdof   <- x$df.residual
aic_n  <- aic/obs
ll     <- xvars - aic/2
bic_r  <- x$deviance - (rdof * log(obs))
bic_l  <- -2*ll + xvars * log(obs)
bic_qh <- -2*(ll - xvars * log(xvars))/obs
return(list("AIC" = aic, "AICn" = aic_n, "BIC" = bic_l, "BICqh" = bic_qh))
}
