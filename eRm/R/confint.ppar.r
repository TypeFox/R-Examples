confint.ppar <- function(object, parm, level = 0.95, ...)
{
#parm...either "beta" or "eta"
#object of class "ppar"

a <- (1 - level)/2
a <- c(a, 1 - a)
pct <- paste(a*100,"%")
fac <- qnorm(a)

cf <- object$thetapar
ses <- object$se.theta

ci <- list(NULL)
for (i in 1:length(cf)) {
  ci[[i]] <- array(NA, dim = c(length(cf[[i]]), 2), dimnames = list(names(object$thetapar[[i]]),pct))
  ci[[i]][] <- cf[[i]] + ses[[i]] %o% fac
}
names(ci) <- paste("NAgroup",1:length(ci),sep="")

ci
}