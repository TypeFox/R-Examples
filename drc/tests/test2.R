## Test provided by Pamela Hutchinson 2007-02-05

#pam <- read.csv("c://stat//projects//R//drcurves//r-part//pkgfolder//drc//tests//test2.redroot_dose.csv")
pam <- read.csv("test2.redroot_dose.csv")

library(drc)

## Initial model: different parameters for the resistant curve
##  and for the susceptible curve

## Temporarily taken out (May 23 2009)

#m1 <- drm(biomass~dose, population, data=pam, fct=LL.4(method = "3"))
#summary(m1)
#plot(m1)


## Reduced model: common slope, lower and upper parameter for both curves,
##  but e/ED50 parameters differ

#m2 <- drm(biomass~dose, population, data=pam, fct=LL.4(method = "3"), pmodels=data.frame(1,1,1,population))
#summary(m2)
#plot(m2)
#anova(m2,m1)

## Further reduced model: all parameters in common
##  This model is too simple: it is rejected

#m3 <- drm(biomass~dose, population, data=pam, fct=LL.4(method = "3"), pmodels=data.frame(1,1,1,1))
#summary(m3)
#plot(m3)
#anova(m3,m2)


## Final model is 'm2' with different ED50 values for resistant and susceptible,
##  but the remaining parameters in common

