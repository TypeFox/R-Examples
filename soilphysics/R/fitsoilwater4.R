fitsoilwater4 <-
function(theta, psi, Bd, model = c("Silva", "Ross"))
{
      model <- match.arg(model)
      if (model == "Silva") {
         fit1. <- lm(log(theta) ~ Bd + log(psi))
         a. <- coef(fit1.)
         fit <- nls(theta ~ exp(a + b*Bd) * psi^c,
            start = list(a = a.[1], b = a.[2], c = a.[3]))
      } else {
         fit1. <- lm(log(theta) ~ log(psi))
         a. <- coef(fit1.)
         fit <- nls(theta ~ a * psi^b,
            start = list(a = exp(a.[1]), b = a.[2]))
      }
      return(fit)
}
