### R code from vignette source 'JSS-gsbDesign.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
remove(list=ls())
set.seed(144)
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
if(!file.exists("./Figures/"))
    dir.create("./Figures/",recursive=TRUE)


###################################################
### code chunk number 2: instloading
###################################################
library("gsbDesign")


###################################################
### code chunk number 3: des1
###################################################
design1 <- gsbDesign(nr.stages = 4, patients = c(10, 20), 
  sigma = c(7, 7), criteria.success = c(0, 0.8, 7, 0.5), 
  criteria.futility = c(2, 0.8), prior.difference = c(3, 5, 2))


###################################################
### code chunk number 4: printDes1
###################################################
names(design1)


###################################################
### code chunk number 5: des2
###################################################
design2 <- gsbDesign(nr.stages = 4, patients = c(10, 20), sigma = c(7, 7),
  criteria.success = c(0, 0.8, 7, 0.5), criteria.futility = c(2, 0.8),
  prior.control = c(3, 5), prior.treatment = c(6, 2))


###################################################
### code chunk number 6: sim1
###################################################
simulation1 <- gsbSimulation(truth = c(-10, 20, 60), 
  type.update = "treatment effect", method = "numerical integration")


###################################################
### code chunk number 7: sim1a
###################################################
simulation1a <- gsbSimulation(truth = c(-10, 20, 60), 
  type.update = "treatment effect", method = "simulation", 
  nr.sim = 50000, warnings.sensitivity = 100, seed = "generate")


###################################################
### code chunk number 8: sim2
###################################################
simulation2 <- gsbSimulation(truth = list(seq(-5, 5, 3), seq(0, 5, 3)),
  type.update = "per arm", method = "simulation", grid.type = "table",
  nr.sim = 10000, warnings.sensitivity = 500, seed = "generate")


###################################################
### code chunk number 9: gsb1
###################################################
oc1 <- gsb(design1, simulation1)


###################################################
### code chunk number 10: gsb1desc
###################################################
names(oc1)


###################################################
### code chunk number 11: gsb1summm
###################################################
summ.oc1 <- summary(oc1, atDelta = c(0, 2, 7))


###################################################
### code chunk number 12: gsb1summ
###################################################
summary(oc1, atDelta = c(0, 2, 7))


###################################################
### code chunk number 13: gsb1fig1prep (eval = FALSE)
###################################################
## plot(oc1, what = "cumulative all")


###################################################
### code chunk number 14: gsb1fig1
###################################################
p1 <- plot(oc1, what = "cumulative all")
print(p1)


###################################################
### code chunk number 15: gsb1fig1values (eval = FALSE)
###################################################
## ## tricks Sweave because of ugly formatting
## formals(gsbDesign:::plot.gsbMainOut)$what


###################################################
### code chunk number 16: gsb1fig2prep (eval = FALSE)
###################################################
## plot(oc1, what = "sample size")
## plot(oc1, what = "boundary")


###################################################
### code chunk number 17: gsb1fig2
###################################################
p2 <- plot(oc1, what = "sample size")
p3 <- plot(oc1, what = "boundary")
print(p2, position = c(0,0, 0.5, 1), more = TRUE)
print(p3, position = c(0.5, 0,1, 1), more = FALSE)


###################################################
### code chunk number 18: gsb1tab1
###################################################
tab(oc1, what = "cumulative success", atDelta = c(0, 2, 7), digits = 4, 
  export = FALSE)


###################################################
### code chunk number 19: gsb2
###################################################
oc2 <- gsb(design2, simulation2)


###################################################
### code chunk number 20: gsb2tab1
###################################################
tab(oc2, what = "sample size", digits = 0)


###################################################
### code chunk number 21: gsb2fig1
###################################################
simulation2a <- gsbSimulation(truth = c(-5, 5, 0, 5, 50), 
  type.update = "per arm", method = "simulation", grid.type = "plot",
  nr.sim = 100000, warnings.sensitivity = 50, seed = "generate")
oc2a <- gsb(design2, simulation2a)
p12 <- plot(oc2a, what = "cumulative all", delta.grid = FALSE)
print(p12)


###################################################
### code chunk number 22: preliminaries
###################################################
remove(list=ls())
set.seed(155)


###################################################
### code chunk number 23: desPoC1
###################################################
desPoC1 <- gsbDesign(nr.stages = 1, patients = c(40, 40), 
  sigma = c(88, 88), criteria.success = c(0, 0.95), 
  criteria.futility = c(NA, NA))
simPoC1 <- gsbSimulation(truth = c(-50, 100, 60),
  type.update = "treatment effect", method = "numerical integration")
ocPoC1 <- gsb(desPoC1, simPoC1)
summary(ocPoC1, atDelta = c(0, 50))


###################################################
### code chunk number 24: desPoC2
###################################################
desPoC2 <- gsbDesign(nr.stages = 1, patients = c(20, 20), 
  sigma = c(88, 88), criteria.success = c(0, 0.975, 50, 0.5), 
  criteria.futility = c(40, 0.9))
simPoC2 <- gsbSimulation(truth = c(0, 70, 60),
  type.update = "treatment effect", method = "numerical integration")
ocPoC2 <- gsb(desPoC2, simPoC2)
summary(ocPoC2, atDelta = c(0, 40, 50, 60))


###################################################
### code chunk number 25: desPoC4
###################################################
desPoC4 <- gsbDesign(nr.stages = 2, patients = c(20, 20), 
  sigma = c(88, 88), criteria.success = c(0, 0.975, 50, 0.5), 
  criteria.futility = c(40, 0.9))
simPoC4 <- gsbSimulation(truth = c(0, 70, 60),
  type.update = "treatment effect", method = "numerical integration")
ocPoC4 <- gsb(desPoC4, simPoC4)
summary(ocPoC4, atDelta = c(0, 40, 50, 60, 70))


###################################################
### code chunk number 26: valuesDesPoC4
###################################################
SPoC4 <- summary(ocPoC4, atDelta = c(0, 40, 50, 60, 70))


###################################################
### code chunk number 27: desPoC5
###################################################
desPoC5 <- gsbDesign(nr.stages = 2, patients = c(10, 20), 
  sigma = c(88, 88), criteria.success = c(0, 0.975, 50, 0.5), 
  criteria.futility = c(40, 0.9), prior.control = c(49, 20))
simPoC5 <- gsbSimulation(truth = cbind(rep(c(30, 50, 70), each = 5),
  c(30, 70, 80, 90, 100, 50, 90, 100, 110, 120, 70, 110, 120, 130, 140)),
  nr.sim = 20000, type.update = "per arm", method = "simulation", 
  grid.type = "manually")
ocPoC5 <- gsb(desPoC5, simPoC5)


###################################################
### code chunk number 28: maketabPoC5
###################################################
tS1 <- subset(ocPoC5$OC,type=="cumulative success" & delta.control==50 & stage=="stage 1")[,c("delta","value")]
tS2 <- subset(ocPoC5$OC,type=="cumulative success" & delta.control==50 & stage=="stage 2")[,"value"]
tF1 <- subset(ocPoC5$OC,type=="cumulative futility" & delta.control==50 & stage=="stage 1")[,"value"]
tF2 <- subset(ocPoC5$OC,type=="cumulative futility" & delta.control==50 & stage=="stage 2")[,"value"]
tN <- ceiling(subset(ocPoC5$OC,type=="sample size" & delta.control==50 & stage=="stage 2")[,"value"])
tSum <- cbind(tS1,tF1,tS2,tF2,tN)
colnames(tSum) <- c("$\\delta$","Interim\nsuccess","Interim\nfutility","Final\nsuccess","Final\nfutility","Expected N")


###################################################
### code chunk number 29: tabPoC5
###################################################
print(xtable(tSum,caption="Operating characteristics of the two stage design.",
 label = "tab1", align=rep("c",7),digits=c(0,0,3,3,3,3,0)),
     stype="latex",sanitize.colnames.function=function(x){x},
     include.rownames=FALSE)


