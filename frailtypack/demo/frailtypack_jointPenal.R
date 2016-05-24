# Virginie Rondeau 2012-4-12 for optimx

options(digits=12)
if(!require("frailtypack"))stop("this test requires package frailtypack.")
if(!require("survival"))stop("this test requires survival.")
if(!require("boot"))stop("this test requires boot.")
if(!require("MASS"))stop("this test requires MASS.")
if(!require("survC1"))stop("this test requires survC1.")

cat("frailtypack test for joint model ...")

#########################################################################
### JOINT frailty model with gap times
#########################################################################

data("readmission")

modJoint.gap <- frailtyPenal(Surv(time, event) ~ cluster(id) +
                               dukes + charlson + sex + chemo + terminal(death),
                             formula.terminalEvent = ~ dukes + charlson + sex + chemo,
                             data = readmission , n.knots = 8 , kappa = c(2.11e+08,9.53e+11))


print(modJoint.gap, digits = 4)

### Print the hazard ratios
summary(modJoint.gap, level = 0.95)

#########################################################################
### Stratified JOINT frailty model with gap times
#########################################################################

modJoint.str <- frailtyPenal(Surv(time, event) ~ cluster(id) +
                               dukes + charlson + strata(sex) + chemo + terminal(death),
                             formula.terminalEvent = ~ dukes + charlson + sex + chemo,
                             data = readmission, n.knots = 8, kappa = c(2.11e+08,2.11e+08,9.53e+11))

print(modJoint.str, digits = 4)

#########################################################################
### JOINT frailty model without alpha parameter (more flexible)
#########################################################################

modJoint.wa <- frailtyPenal(Surv(time, event) ~ cluster(id) +
                              dukes + charlson + sex + chemo + terminal(death),
                            formula.terminalEvent = ~ dukes + charlson + sex + chemo,
                            data = readmission, n.knots = 8, kappa = c(2.11e+08,9.53e+11), Alpha = "none")

print(modJoint.wa, digits = 4)

#########################################################################
### JOINT frailty model for clustered data
#########################################################################

readmission <- transform(readmission,group=id%%31+1)

modJoint.clus <- frailtyPenal(Surv(t.start, t.stop, event) ~ cluster(group) + num.id(id) +
                                dukes + charlson + sex + chemo + terminal(death),
                              formula.terminalEvent = ~ dukes + charlson + sex + chemo,
                              data = readmission, recurrentAG = TRUE,  n.knots = 10, kappa = c(2.11e+08,9.53e+11))

print(modJoint.clus, digits = 4)

########################################################################
### Figures
########################################################################
pdf(file="fig3.pdf", height = 3.6, width = 8.1)
par(mfrow = c(1, 3))
plot(modJoint.gap, type.plot = "survival", event = "recurrent", main = "Recurrent",
     conf.bands = TRUE, pos.legend = "topleft", cex.legend = 1.2, ylim = c(0, 1.2))
plot(modJoint.gap, type.plot = "survival", event = "terminal", main = "Terminal",
     conf.bands = TRUE, pos.legend = "topleft", cex.legend = 1.2, ylim = c(0, 1.2))
plot(modJoint.gap, type.plot = "survival", event = "both", main = "Both",
     conf.bands = TRUE, pos.legend = "topleft", cex.legend = 1, ylim = c(0, 1.2))



################################################################################################################
###  General Joint model (recurrent and terminal events) with 2 covariates
################################################################################################################

library("frailtypack")
data(readmission)
modJoint.general <- frailtyPenal(Surv(time,event) ~ cluster(id) + dukes +
                             charlson + sex  + chemo + terminal(death),formula.terminalEvent = ~ dukes +
                             charlson + sex + chemo, data = readmission, jointGeneral = TRUE,  n.knots = 8,
                             kappa = c(2.11e+08, 9.53e+11))

# Start writing to frailtyPack_Jgeneral_Results.txt file
sink("frailtyPack_Jgeneral_Results.txt")
print(modJoint.general, digits = 4)
# Stop writing to the file
sink()
