library(stocc)
data(habData)
data(visitData)

names <- list(
  visit=list(site="site",obs="obs"),
  site=list(site="site", coords=c("x","y"))
  )

## This takes a while...
fit <- spatial.occupancy(
  detection.model = ~ detCov1 + detCov2,
  occupancy.model = ~ habCov1 + habCov2,
  spatial.model = list(model="rsr", threshold=1, moran.cut=0.1*nrow(habData)), 
  so.data = make.so.data(visitData, habData, names),
  prior = list(a.tau=0.5, b.tau=0.00005, Q.b=0.1, Q.g=0.1),
  control = list(burnin=1000/5, iter=4000/5, thin=5)
  )

## Model fit -- Minimum posterior predictive loss
fit$D.m # = G.m + P.m

fit$G.m # Goodness of fit
fit$P.m # Complexity penalty


## Plot of detection parameter MCMC
plot(fit$beta)

## 95\% CI for detection parameters
HPDinterval(fit$beta, 0.9)

## Plot of occupancy parameters.
dev.new()
plot(fit$gamma)

## 95\% CI for occupancy parameters
HPDinterval(fit$gamma, 0.9)


## Plot of spatial variance parameter.
dev.new()
plot(fit$tau)

## Spatial plot of the occupancy process. Circles represent sites with observed occupancy and dots represent sites for which
## occupancy was never confirmed
dev.new()
attach(fit$occupancy.df)
image(x=seq(0.5,39.5,1), y=seq(0.5,39.5,1), z=t(matrix(psi.est, 40)), main="Estimated occupancy process", xlab="x", ylab="y")
points(x[samp==1], y[samp==1], cex=0.25, col="blue", pch=20)
points(x[real.occ.est==1 & samp==1], y[real.occ.est==1 & samp==1], cex=1, col="blue", pch=1)
detach(fit$occupancy.df)

## Plot of the spatial random effect process. Circles represent sites with observed occupancy and dots represent sites for which
## occupancy was never confirmed
dev.new()
attach(fit$occupancy.df)
image(x=seq(0.5,39.5,1), y=seq(0.5,39.5,1), z=t(matrix(eta.est, 40)), main="Estimated spatial effect process", xlab="x", ylab="y")
points(x[samp==1], y[samp==1], cex=0.25, col="blue", pch=20)
points(x[real.occ.est==1 & samp==1], y[real.occ.est==1 & samp==1], cex=1, col="blue", pch=1)
detach(fit$occupancy.df)

