### R code from vignette source 'examples.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(prompt = "R> ")
options(SweaveHooks = list(
  cex = function() par(cex.lab = 1.3, cex.axis = 1.3)))


###################################################
### code chunk number 2: FOCUS_2006_L1_data
###################################################
library("kinfit")
FOCUS_2006_L1 = kinobject("Parent", "Degradation data", "")
FOCUS_2006_L1$data = data.frame(
  t = rep(c(0, 1, 2, 3, 5, 7, 14, 21, 30), each = 2),
  parent = c(88.3, 91.4, 85.6, 84.5, 78.9, 77.6, 
             72.0, 71.9, 50.3, 59.4, 47.0, 45.1,
             27.7, 27.3, 10.0, 10.4, 2.9, 4.0))


###################################################
### code chunk number 3: L1
###################################################
FOCUS_2006_L1$fits <- kinfit(FOCUS_2006_L1$data, 
  kinmodels = c("SFO", "FOMC", "DFOP"))
FOCUS_2006_L1$results <- kinresults(FOCUS_2006_L1$fits)
kinreport(FOCUS_2006_L1)


###################################################
### code chunk number 4: L1_2
###################################################
FOCUS_2006_L1$fits <- kinfit(FOCUS_2006_L1$data, 
  kinmodels = c("SFO", "FOMC", "DFOP"),
  start.FOMC = list(parent.0 = 92.47, alpha = 1.35e11, beta = 1.41e12))
FOCUS_2006_L1$results <- kinresults(FOCUS_2006_L1$fits)
kinreport(FOCUS_2006_L1)


###################################################
### code chunk number 5: L1_SFO_plot
###################################################
kinplot(FOCUS_2006_L1, ylab = "Observed")


###################################################
### code chunk number 6: L1_SFO_residuals
###################################################
kinresplot(FOCUS_2006_L1, "SFO", ylab = "Observed")


###################################################
### code chunk number 7: FOCUS_2006_L2_data
###################################################
FOCUS_2006_L2 = kinobject("Parent", "Degradation data", "")
FOCUS_2006_L2$data = data.frame(
  t = rep(c(0, 1, 3, 7, 14, 28), each = 2),
  parent = c(96.1, 91.8, 41.4, 38.7,
             19.3, 22.3, 4.6, 4.6,
             2.6, 1.2, 0.3, 0.6))


###################################################
### code chunk number 8: L2
###################################################
FOCUS_2006_L2$fits <- kinfit(FOCUS_2006_L2$data, 
  kinmodels = c("SFO", "FOMC", "DFOP"))
FOCUS_2006_L2$results <- kinresults(FOCUS_2006_L2$fits)
kinreport(FOCUS_2006_L2)


###################################################
### code chunk number 9: L2_2
###################################################
FOCUS_2006_L2$fits <- kinfit(FOCUS_2006_L2$data, 
  kinmodels = c("SFO", "FOMC", "DFOP"),
  start.DFOP = list(parent.0 = 94, g = 0.4, k1 = 142, k2 = 0.34))
FOCUS_2006_L2$results <- kinresults(FOCUS_2006_L2$fits)
kinreport(FOCUS_2006_L2)


###################################################
### code chunk number 10: L2_plot
###################################################
kinplot(FOCUS_2006_L2, ylab = "Observed")


###################################################
### code chunk number 11: L2_resplot
###################################################
par(mfrow=c(2,1))
kinresplot(FOCUS_2006_L2, "SFO", ylab = "Observed")
kinresplot(FOCUS_2006_L2, "FOMC", ylab = "Observed")


###################################################
### code chunk number 12: FOCUS_2006_L3
###################################################
FOCUS_2006_L3 = kinobject("Parent", "Degradation data", "")
FOCUS_2006_L3$data = data.frame(
  t = c(0, 3, 7, 14, 30, 60, 91, 120),
  parent = c(97.8, 60, 51, 43, 35, 22, 15, 12))
FOCUS_2006_L3$fits <- kinfit(FOCUS_2006_L3$data, 
  kinmodels = c("SFO", "FOMC", "DFOP"))
FOCUS_2006_L3$results <- kinresults(FOCUS_2006_L3$fits)
kinreport(FOCUS_2006_L3)


###################################################
### code chunk number 13: FOCUS_2006_L3_2
###################################################
FOCUS_2006_L3$fits <- kinfit(FOCUS_2006_L3$data, 
  kinmodels = c("SFO", "FOMC", "DFOP"),
  start.FOMC = list(parent.0 = 100, alpha = 0.5, beta = 2))
FOCUS_2006_L3$results <- kinresults(FOCUS_2006_L3$fits)
kinreport(FOCUS_2006_L3)
kinplot(FOCUS_2006_L3, ylab = "Observed")


###################################################
### code chunk number 14: FOCUS_2006_L4
###################################################
FOCUS_2006_L4 = kinobject("Parent", "Degradation data", "")
FOCUS_2006_L4$data = data.frame(
  t = c(0, 3, 7, 14, 30, 60, 91, 120),
  parent = c(96.6, 96.3, 94.3, 88.8, 74.9, 59.9, 53.5, 49.0))
FOCUS_2006_L4$fits <- kinfit(FOCUS_2006_L4$data, 
  kinmodels = c("SFO", "FOMC", "DFOP"))
FOCUS_2006_L4$results <- kinresults(FOCUS_2006_L4$fits)
kinreport(FOCUS_2006_L4)
kinplot(FOCUS_2006_L4, ylab = "Observed")


