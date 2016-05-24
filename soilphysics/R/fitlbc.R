fitlbc <-
function(theta, sigmaP, ...)
{
   if (length(theta) != length(sigmaP))
      stop ("incompatible dimensions!")
   dat <- data.frame(theta, sigmaP)
   ini <- try( getInitiallbc(theta, sigmaP) )
   if (inherits(ini, "try-error")) {
      ini <- coef(lm(log10(sigmaP) ~ theta, data = dat))
   }
   fit <- nls(sigmaP ~ SSlbc(theta, b0, b1), data = dat,
      start = list(b0 = ini[1], b1 = ini[2]), ...)
   return(fit)
}
