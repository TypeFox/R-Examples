fitbusscher <-
function(Pr, theta, Bd, ...)
{
   if (length(Pr) != length(theta) ||
   length(theta) != length(Bd))
      stop ("incompatible dimensions!")
   dat <- data.frame(Pr, theta, Bd)

   ini <- lm(log10(Pr) ~ log10(theta) + log10(Bd),
      data = dat)$coef

   fit <- nls(Pr ~ b0 * (theta ^ b1) * (Bd ^ b2),
      start = list(b0 = 10^ini[1], b1 = ini[2], b2 = ini[3]),
      data = dat, ...)
   return(fit)
}
