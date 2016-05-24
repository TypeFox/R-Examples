# Virginie Rondeau 2012-4-12 for optimx

options(digits=12)
if(!require("frailtypack"))stop("this test requires package frailtypack.")
if(!require("survival"))stop("this test requires survival.")
if(!require("boot"))stop("this test requires boot.")
if(!require("MASS"))stop("this test requires MASS.")
if(!require("survC1"))stop("this test requires survC1.")

cat("frailtypack test for shared model ...")


########################################################################
### COX proportionnal hazard model with gap times
########################################################################

data("readmission")

mod.cox.gap <- frailtyPenal(Surv(time, event) ~ dukes +
  charlson + sex + chemo, n.knots = 10, kappa = 1,
  data = readmission, cross.validation = TRUE)

print(mod.cox.gap, digits = 4)

########################################################################
### Shared frailty model with gap times
########################################################################

mod.sha.gap <- frailtyPenal(Surv(time, event) ~ cluster(id) +
  dukes + charlson + sex + chemo, n.knots = 10, kappa = 1,
  data = readmission, cross.validation = TRUE)

print(mod.sha.gap, digits = 4)

########################################################################
### Shared frailty model with log-normal distribution of frailties
########################################################################

mod.sha.log <- frailtyPenal(Surv(time, event) ~ cluster(id) +
  dukes + charlson + sex + chemo, n.knots = 10, kappa = 1,
  data = readmission, cross.validation = TRUE, RandDist = "LogN")

print(mod.sha.log, digits = 4)

#########################################################################
### Stratified shared frailty model with gap times
#########################################################################

mod.sha.str.gap <- frailtyPenal(Surv(time, event) ~ cluster(id) +
  charlson + dukes + chemo + strata(sex), n.knots = 10,
  kappa = c(2.11e+08,2.11e+08), data = readmission)

print(mod.sha.str.gap, digits = 4)

#########################################################################
### Shared frailty model with time-varying effect of covariates
#########################################################################

mod.sha.time <- frailtyPenal(Surv(time, event) ~ cluster(id) +
  dukes + charlson + timedep(sex) + chemo, n.knots = 8,
  kappa = 1, data = readmission, betaknots = 3, betaorder = 1)

print(mod.sha.time, digits = 4)

#########################################################################
### Shared frailty model with interval-censored data
#########################################################################

data(bcos)
bcos$event <- ifelse(bcos$left!=bcos$right,1,0)
bcos$group <- c(rep(1:20,4),1:14)

mod.sha.ic <- frailtyPenal(SurvIC(left, right, event) ~ cluster(group) +
  treatment, n.knots = 8, kappa = 10000, data = bcos)

print(mod.sha.ic, digits = 4)

########################################################################
### Figures
########################################################################

pdf(file = "fig1.pdf")
par(mfrow = c(1, 3))
plot(mod.cox.gap, type.plot = "survival", main = "Cox model", conf.bands = TRUE)
plot(mod.sha.gap, type.plot = "survival", main = "Shared", conf.bands = TRUE)
plot(mod.sha.str.gap, type.plot = "survival", main = "Shared + Stratification",
  conf.bands = TRUE, pos.legend = "bottomleft", cex.legend = 1)
