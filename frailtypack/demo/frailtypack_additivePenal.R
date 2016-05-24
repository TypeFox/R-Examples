# Virginie Rondeau 2012-4-12 for optimx

options(digits=12)
if(!require("frailtypack"))stop("this test requires package frailtypack.")
if(!require("survival"))stop("this test requires survival.")
if(!require("boot"))stop("this test requires boot.")
if(!require("MASS"))stop("this test requires MASS.")
if(!require("survC1"))stop("this test requires survC1.")

cat("frailtypack test for additive model ...")


########################################################################
### ADDITIVE frailty model with no correlation between random effects
########################################################################

data("dataAdditive")

modAdd2cov.withoutCorr <- additivePenal(Surv(t1, t2, event) ~ cluster(group) +
  var1 + var2 + slope(var1), cross.validation = TRUE,
  correlation = FALSE, data = dataAdditive, n.knots = 10, kappa = 1)

print(modAdd2cov.withoutCorr, digits = 4)
plot(modAdd2cov.withoutCorr)

########################################################################
### ADDITIVE frailty model with a correlation between random effects
########################################################################

modAdd2cov.withCorr <- additivePenal(Surv(t1,t2,event) ~ cluster(group) +
  var1 + var2 + slope(var1), cross.validation = TRUE,
  data = dataAdditive, correlation = TRUE, n.knots = 10, kappa = 1)

print(modAdd2cov.withCorr, digits = 4)
plot(modAdd2cov.withCorr)


########################################################################
### Figures
########################################################################

pdf(file = "fig2.pdf", height = 5.5, width = 9)
par(mfrow = c(1,2))
plot(modAdd2cov.withoutCorr, type.plot = "hazard", main = "Correlation=False",
  conf.bands = TRUE)
plot(modAdd2cov.withCorr, type.plot = "hazard", main = "Correlation=True",
  conf.bands = TRUE)