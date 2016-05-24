## Case Study: Design of a Proof-of-Concept (PoC) trial in Crohn's disease
## See Section 5 of the vignette JSS-gsbDesign.pdf for more information.

set.seed(155)
desPoC1 <- gsbDesign(nr.stages = 1, patients = c(40, 40), 
  sigma = c(88, 88), criteria.success = c(0, 0.95), 
  criteria.futility = c(NA, NA))
simPoC1 <- gsbSimulation(truth = c(-50, 100, 60),
  type.update = "treatment effect", method = "numerical integration")
ocPoC1 <- gsb(desPoC1, simPoC1)
summary(ocPoC1, atDelta = c(0, 50))

desPoC2 <- gsbDesign(nr.stages = 1, patients = c(20, 20), 
  sigma = c(88, 88), criteria.success = c(0, 0.975, 50, 0.5), 
  criteria.futility = c(40, 0.9))
simPoC2 <- gsbSimulation(truth = c(0, 70, 60),
  type.update = "treatment effect", method = "numerical integration")
ocPoC2 <- gsb(desPoC2, simPoC2)
summary(ocPoC2, atDelta = c(0, 40, 50, 60))


desPoC4 <- gsbDesign(nr.stages = 2, patients = c(20, 20), 
  sigma = c(88, 88), criteria.success = c(0, 0.975, 50, 0.5), 
  criteria.futility = c(40, 0.9))
simPoC4 <- gsbSimulation(truth = c(0, 70, 60),
  type.update = "treatment effect", method = "numerical integration")
ocPoC4 <- gsb(desPoC4, simPoC4)
summary(ocPoC4, atDelta = c(0, 40, 50, 60, 70))


SPoC4 <- summary(ocPoC4, atDelta = c(0, 40, 50, 60, 70))


desPoC5 <- gsbDesign(nr.stages = 2, patients = c(10, 20), 
  sigma = c(88, 88), criteria.success = c(0, 0.975, 50, 0.5), 
  criteria.futility = c(40, 0.9), prior.control = c(49, 20))
simPoC5 <- gsbSimulation(truth = cbind(rep(c(30, 50, 70), each = 5),
  c(30, 70, 80, 90, 100, 50, 90, 100, 110, 120, 70, 110, 120, 130, 140)),
  nr.sim = 20000, type.update = "per arm", method = "simulation", 
  grid.type = "manually")
ocPoC5 <- gsb(desPoC5, simPoC5)


## create Table 1
tS1 <- subset(ocPoC5$OC,type=="cumulative success" & delta.control==50 & stage=="stage 1")[,c("delta","value")]
tS2 <- subset(ocPoC5$OC,type=="cumulative success" & delta.control==50 & stage=="stage 2")[,"value"]
tF1 <- subset(ocPoC5$OC,type=="cumulative futility" & delta.control==50 & stage=="stage 1")[,"value"]
tF2 <- subset(ocPoC5$OC,type=="cumulative futility" & delta.control==50 & stage=="stage 2")[,"value"]
tN <- ceiling(subset(ocPoC5$OC,type=="sample size" & delta.control==50 & stage=="stage 2")[,"value"])
tSum <- cbind(tS1,tF1,tS2,tF2,tN)
colnames(tSum) <- c("$\\delta$","Interim\nsuccess","Interim\nfutility","Final\nsuccess","Final\nfutility","Expected N")


