## ---- echo=FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE)

## ------------------------------------------------------------------------
library("bayesmeta")
data("Rubin1981")
print(Rubin1981)

## ------------------------------------------------------------------------
data("CrinsEtAl2014")
print(CrinsEtAl2014[,c(1,10,12,13,15)])

## ---- message=FALSE------------------------------------------------------
library("metafor")
crins.es <- escalc(measure="OR",
                   ai=exp.AR.events,  n1i=exp.total,
                   ci=cont.AR.events, n2i=cont.total,
                   slab=publication, data=CrinsEtAl2014)
print(crins.es[,c("publication", "yi", "vi")])

## ------------------------------------------------------------------------
data("Rubin1981")
taupriordensity <- function(t){dhalfcauchy(t, scale=25)}
schools.example.1 <- bayesmeta(y     = Rubin1981[,"effect"],
                               sigma = Rubin1981[,"stderr"],
                               label = Rubin1981[,"school"],
                               mu.prior.mean=0, mu.prior.sd=50,
                               tau.prior=taupriordensity)

## ------------------------------------------------------------------------
print(schools.example.1)

## ---- eval=FALSE---------------------------------------------------------
#  plot(schools.example.1, prior=TRUE)

## ---- fig.width=6.0, fig.height=7.0, echo=FALSE--------------------------
par(mfrow=c(2,2))
plot(schools.example.1, prior=TRUE)
par(mfrow=c(1,1))

## ------------------------------------------------------------------------
schools.example.2 <- bayesmeta(y     = Rubin1981[,"effect"],
                               sigma = Rubin1981[,"stderr"],
                               label = Rubin1981[,"school"])

## ------------------------------------------------------------------------
print(schools.example.1$summary)
print(schools.example.2$summary)

## ---- fig.width=5.0, fig.height=5.0--------------------------------------
# evaluate posterior densities:
x <- seq(from=-10, to=30, length=100)
plot(x, schools.example.1$dposterior(mu=x), type="l", col="red",
     xlab=expression("effect "*mu), ylab="posterior density")
lines(x, schools.example.2$dposterior(mu=x), type="l", col="blue", lty="dashed")
abline(h=0, col="darkgrey")

## ------------------------------------------------------------------------
# posterior probability of mu > 0:
1 - schools.example.1$pposterior(mu=0)

## ------------------------------------------------------------------------
# 95% posterior upper limit on the effect mu:
schools.example.1$qposterior(mu.p=0.95)

## ------------------------------------------------------------------------
# 95% posterior upper limit on the heterogeneity tau:
schools.example.1$qposterior(tau.p=0.95)

## ------------------------------------------------------------------------
# 95% credible intervals for the effect mu:
schools.example.1$post.interval(mu.level=0.95)

## ------------------------------------------------------------------------
# 95% credible intervals for the effect mu:
schools.example.1$post.interval(tau.level=0.95)
schools.example.1$post.interval(tau.level=0.95, method="central")
schools.example.1$qposterior(tau.p=c(0.025, 0.975))

