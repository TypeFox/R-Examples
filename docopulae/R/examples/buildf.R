## for an actual use case see examples for param

library(copula)
library(mvtnorm)

## build bivariate normal pdf
margins = function(y, theta) {
    mu = c(theta$mu1, theta$mu2)
    cbind(dnorm(y, mu), pnorm(y, mu))
}
copula = copula::normalCopula()

f = buildf(margins, copula, 'alpha')
f

## plot density
theta = list(mu1=2, mu2=-3, alpha=0.4)
y1 = seq(0, 4, length.out=51)
y2 = seq(-5, -1, length.out=51)
z = outer(y1, y2, function(y1, y2) apply(cbind(y1, y2), 1, f, theta))
contour(y1, y2, z)

## add theoretical density
copula@parameters = theta$alpha
z2 = outer(y1, y2, function(y1, y2)
    dmvnorm(cbind(y1 - theta$mu1, y2 - theta$mu2), sigma=getSigma(copula)))
contour(y1, y2, z2, col='red', add=TRUE)

## build bivariate pdf with normal margins and clayton copula
margins = list(alist(pdf=dnorm(y1, mu1, 1),
                     cdf=pnorm(y1, mu1, 1)),
               alist(pdf=dnorm(y2, mu2, 1),
                     cdf=pnorm(y2, mu2, 1)))
copula = claytonCopula()
ff = buildf(margins, copula)
f = expr2f(ff, yMap=list(y1=1, y2=2),
               thetaMap=list(mu1='mu1', mu2='mu2', alpha='alpha'))
f

margins = function(y, theta) {
    mu = c(theta$mu1, theta$mu2)
    cbind(dnorm(y, mu, 1), pnorm(y, mu, 1))
}
f2 = buildf(margins, copula, 'alpha')
f2

## plot both densities
theta = list(mu1=2, mu2=-3, alpha=2) # tau = 0.5

y1 = seq(0, 4, length.out=51)
y2 = seq(-5, -1, length.out=51)
z = outer(y1, y2, function(y1, y2) apply(cbind(y1, y2), 1, f, theta))
contour(y1, y2, z)

z2 = outer(y1, y2, function(y1, y2) apply(cbind(y1, y2), 1, f2, theta))
contour(y1, y2, z2, col='red', add=TRUE)
