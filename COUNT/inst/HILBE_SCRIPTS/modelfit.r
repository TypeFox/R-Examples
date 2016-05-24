# modelfit function to calc AIC and BIC statistics post estimation
# Joseph M. Hilbe 12January, 2010
modelfit  <- function(x)  {
obs    <- x$df.null + 1
aic    <- x$aic
xvars  <- x$rank
rdof   <- x$df.residual
aic_n  <- aic/obs
ll     <- xvars - aic/2
bic_r  <- x$deviance - (rdof * log(obs))
bic_l  <- -2*ll + xvars * log(obs)
bic_qh <- -2*(ll - xvars * log(xvars))/obs
return(list("AIC" = aic, "AICn" = aic_n, "BIC" = bic_l, "BICqh" = bic_qh))
}
modelfit(x)  # substitute fitted model name for x
