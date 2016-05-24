confint.eRm <- function(object, parm="beta", level = 0.95, ...)
{
#parm...either "beta" or "eta"
#object of class "eRm"

a <- (1 - level)/2
a <- c(a, 1 - a)
pct <- paste(a*100,"%")
fac <- qnorm(a)

if (parm=="eta") {
  cf <- object$etapar
  ses <- object$se.eta
  dn <- names(object$etapar) }
if (parm=="beta") {
  cf <- object$betapar
  ses <- object$se.beta
  dn <- names(object$betapar) }

ci <- array(NA, dim = c(length(cf), 2), dimnames = list(dn,pct))
ci[] <- cf + ses %o% fac
ci

}



