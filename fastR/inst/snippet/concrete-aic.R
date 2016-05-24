# these two methods give different numerical values
AIC(concrete.lm0)
AIC(concrete.lm1)
extractAIC(concrete.lm0)
extractAIC(concrete.lm1)
# and neither agrees with our definition
aic0 <- 2 * 3 + 9 * log(sum(resid(concrete.lm0)^2)); aic0
aic1 <- 2 * 2 + 9 * log(sum(resid(concrete.lm1)^2)); aic1
# but differences between models are equivalent
aic0 - aic1
AIC(concrete.lm0) - AIC(concrete.lm1)
extractAIC(concrete.lm0)[2] - extractAIC(concrete.lm1)[2]
