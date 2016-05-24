confint.threshold <- function(object, parm, level = 0.95, ...)
{
#object of class "threshold"

a <- (1 - level)/2
a <- c(a, 1 - a)
pct <- paste(a*100,"%")
fac <- qnorm(a)

cf <- object$threshpar
ses <- object$se.thresh
dn <- names(object$threshpar)

ci <- array(NA, dim = c(length(cf), 2), dimnames = list(dn,pct))
ci[] <- cf + ses %o% fac
ci

}