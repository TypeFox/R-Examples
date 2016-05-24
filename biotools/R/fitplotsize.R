fitplotsize <-
function(plotsize, CV)
{
   if (length(plotsize) != length(CV))
      stop ("incompatible dimensions!")

   ini <- lm(log(CV) ~ log(plotsize))$coef

   fit <- nls(CV ~ a * plotsize^(-b),
      start = list(a = ini[1], b = ini[2]))
   print(summary(fit))
   invisible(fit)
}
