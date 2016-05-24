pbsize <- function (kp, gamma=4.5, p=0.15, alpha=5e-8, beta=0.2)
# population-based sample size
# alpha=5e-8, beta=0.8, alpha would give 5% genome-wide significance level
# x2alpha = 29.72 (Q=29.7168)
#
# lambda is the NCP from the marginal table
# pi is the pr(Affected|aa)
{
  z1alpha <- qnorm(1-alpha/2)
  z1beta <- qnorm(1-beta)
  q <- 1-p
  pi <- kp/(gamma*p+q)^2
  lambda <- pi*p*q*(gamma-1)^2/(1-kp)
  n <- (z1alpha+z1beta)^2/lambda
  n
}
