###################################################
### code chunk number 38: Case-A-unit-test
###################################################
test.caseA <- function() {
trialParameters <- list(minimumNumberOfEvents = 20,
                        minimumIncreaseInV = 0.2,
                        numberRecruitedEachYear = c(80, 120, 160, 160),
                        followupTime = 3,
                        adminCensoringTime = 7,
                        interimLookTime = c(1, 2, 3, 5, 7),
                        type1ErrorForResponse = 0.05,
                        type2ErrorForResponse = 0.01,
                        glrBoundarySidedness = "one", # one sided or two-sided
                        type1Error = 0.05,
                        type2Error = 0.10,
                        epsType1 = 1/3,
                        epsType2 = 1/3)
trueParameters <- list(p0 = 0.3,
                       p1 = 0.3,
                       pdiffHyp=0.3,
                       theta = list(
                         alpha = 0,
                         beta = 0,
                         gamma = 0),
                       baselineLambda = 0.35,
                       etaHyp = 0.25)
rngSeed <- 9872831
sp23Design <- generateSP23Design(trueParameters, trialParameters)
print(sp23Design)
trialHistory <- exploreSP23Design(sp23Design, numberOfSimulations=25, rngSeed=rngSeed)
result <- analyzeSP23Design(sp23Design, trialHistory)$designSummary
  checkEquals(result[["numberOfTimesH0RIsRejectedAtFirstLook"]], 0)
  checkEquals(result[["numberOfTimesH0RIsRejected"]], 0)
  checkEquals(result[["numberOfTimesStoppedForFutility"]], 25)
  checkEquals(result[["numberOfTimesH0SIsAccepted"]], 0)
  checkEquals(result[["numberOfTimesH0SIsRejected"]], 0)
  checkEquals(result[["numberOfTimesFutilityDecidedAtLastLook"]], 0)
  checkEquals(result[["numberOfTimesTrialEndedAtLook"]],
          {out <- c(9, 13, 3, 0, 0); names(out) <- trialParameters$interimLookTime; out})

}


###################################################
### code chunk number 39: Case-B-unit-test
###################################################
test.caseB <- function() {
trialParameters <- list(minimumNumberOfEvents = 20,
                        minimumIncreaseInV = 0.2,
                        numberRecruitedEachYear = c(80, 120, 160, 160),
                        followupTime = 3,
                        adminCensoringTime = 7,
                        interimLookTime = c(1, 2, 3, 5, 7),
                        type1ErrorForResponse = 0.05,
                        type2ErrorForResponse = 0.01,
                        glrBoundarySidedness = "one", # one sided or two-sided
                        type1Error = 0.05,
                        type2Error = 0.10,
                        epsType1 = 1/3,
                        epsType2 = 1/3)
trueParameters <- list(p0 = 0.3,
                       p1 = 0.3,
                       pdiffHyp=0.3,
                       theta = list(
                         alpha = log(0.5),
                         beta = 0,
                         gamma = 0),
                       baselineLambda = 0.35,
                       etaHyp = 0.25)
rngSeed <- 324612
sp23Design <- generateSP23Design(trueParameters, trialParameters)
print(sp23Design)
trialHistory <- exploreSP23Design(sp23Design, numberOfSimulations=25, rngSeed=rngSeed)
result <- analyzeSP23Design(sp23Design, trialHistory)$designSummary
  checkEquals(result[["numberOfTimesH0RIsRejectedAtFirstLook"]], 0)
  checkEquals(result[["numberOfTimesH0RIsRejected"]], 0)
  checkEquals(result[["numberOfTimesStoppedForFutility"]], 25)
  checkEquals(result[["numberOfTimesH0SIsAccepted"]], 0)
  checkEquals(result[["numberOfTimesH0SIsRejected"]], 0)
  checkEquals(result[["numberOfTimesFutilityDecidedAtLastLook"]], 0)
  checkEquals(result[["numberOfTimesTrialEndedAtLook"]],
          {out <- c(9, 16, 0, 0, 0); names(out) <- trialParameters$interimLookTime; out})
}


###################################################
### code chunk number 40: Case-C-unit-test
###################################################
test.caseC <- function() {
trialParameters <- list(minimumNumberOfEvents = 20,
                        minimumIncreaseInV = 0.2,
                        numberRecruitedEachYear = c(80, 120, 160, 160),
                        followupTime = 3,
                        adminCensoringTime = 7,
                        interimLookTime = c(1, 2, 3, 5, 7),
                        type1ErrorForResponse = 0.05,
                        type2ErrorForResponse = 0.01,
                        glrBoundarySidedness = "one", # one sided or two-sided
                        type1Error = 0.05,
                        type2Error = 0.10,
                        epsType1 = 1/3,
                        epsType2 = 1/3)
trueParameters <- list(p0 = 0.3,
                       p1 = 0.6,
                       pdiffHyp=0.3,
                       theta = list(
                         alpha = 0,
                         beta = 0,
                         gamma = 0),
                       baselineLambda = 0.35,
                       etaHyp = 0.25)
rngSeed <- 2387487
sp23Design <- generateSP23Design(trueParameters, trialParameters)
print(sp23Design)
trialHistory <- exploreSP23Design(sp23Design, numberOfSimulations=25, rngSeed=rngSeed)
result <- analyzeSP23Design(sp23Design, trialHistory)$designSummary
  checkEquals(result[["numberOfTimesH0RIsRejectedAtFirstLook"]], 14)
  checkEquals(result[["numberOfTimesH0RIsRejected"]], 25)
  checkEquals(result[["numberOfTimesStoppedForFutility"]], 23)
  checkEquals(result[["numberOfTimesH0SIsAccepted"]], 23)
  checkEquals(result[["numberOfTimesH0SIsRejected"]], 2)
  checkEquals(result[["numberOfTimesFutilityDecidedAtLastLook"]], 7)
  checkEquals(result[["numberOfTimesTrialEndedAtLook"]],
          {out <- c(1, 1, 3, 12, 8); names(out) <- trialParameters$interimLookTime; out})
}


###################################################
### code chunk number 41: Case-D-unit-test
###################################################
test.caseD <- function() {
trialParameters <- list(minimumNumberOfEvents = 20,
                        minimumIncreaseInV = 0.2,
                        numberRecruitedEachYear = c(80, 120, 160, 160),
                        followupTime = 3,
                        adminCensoringTime = 7,
                        interimLookTime = c(1, 2, 3, 5, 7),
                        type1ErrorForResponse = 0.05,
                        type2ErrorForResponse = 0.01,
                        glrBoundarySidedness = "one", # one sided or two-sided
                        type1Error = 0.05,
                        type2Error = 0.10,
                        epsType1 = 1/3,
                        epsType2 = 1/3)
trueParameters <- list(p0 = 0.3,
                       p1 = 0.6,
                       pdiffHyp=0.3,
                       theta = list(
                         alpha = log(0.2),
                         beta = log(1.56),
                         gamma = log(1.0)),
                       baselineLambda = 0.35,
                       etaHyp = 0.25)
rngSeed <- 283119
sp23Design <- generateSP23Design(trueParameters, trialParameters)
print(sp23Design)
trialHistory <- exploreSP23Design(sp23Design, numberOfSimulations=25, rngSeed=rngSeed)
result <- analyzeSP23Design(sp23Design, trialHistory)$designSummary
  checkEquals(result[["numberOfTimesH0RIsRejectedAtFirstLook"]], 22)
  checkEquals(result[["numberOfTimesH0RIsRejected"]], 25)
  checkEquals(result[["numberOfTimesStoppedForFutility"]], 24)
  checkEquals(result[["numberOfTimesH0SIsAccepted"]], 24)
  checkEquals(result[["numberOfTimesH0SIsRejected"]], 1)
  checkEquals(result[["numberOfTimesFutilityDecidedAtLastLook"]], 1)
  checkEquals(result[["numberOfTimesTrialEndedAtLook"]],
          {out <- c(0, 3, 10, 11, 1); names(out) <- trialParameters$interimLookTime; out})
}


###################################################
### code chunk number 42: Case-E-unit-test
###################################################
test.caseE <- function() {
trialParameters <- list(minimumNumberOfEvents = 20,
                        minimumIncreaseInV = 0.2,
                        numberRecruitedEachYear = c(80, 120, 160, 160),
                        followupTime = 3,
                        adminCensoringTime = 7,
                        interimLookTime = c(1, 2, 3, 5, 7),
                        type1ErrorForResponse = 0.05,
                        type2ErrorForResponse = 0.01,
                        glrBoundarySidedness = "one", # one sided or two-sided
                        type1Error = 0.05,
                        type2Error = 0.10,
                        epsType1 = 1/3,
                        epsType2 = 1/3)
trueParameters <- list(p0 = 0.3,
                       p1 = 0.6,
                       pdiffHyp=0.3,
                       theta = list(
                         alpha = log(0.2),
                         beta = log(1.73),
                         gamma = log(0.8)),
                       baselineLambda = 0.35,
                       etaHyp = 0.25)
rngSeed <- 98872361
sp23Design <- generateSP23Design(trueParameters, trialParameters)
print(sp23Design)
trialHistory <- exploreSP23Design(sp23Design, numberOfSimulations=25, rngSeed=rngSeed)
result <- analyzeSP23Design(sp23Design, trialHistory)$designSummary
  checkEquals(result[["numberOfTimesH0RIsRejectedAtFirstLook"]], 16)
  checkEquals(result[["numberOfTimesH0RIsRejected"]], 25)
  checkEquals(result[["numberOfTimesStoppedForFutility"]], 25)
  checkEquals(result[["numberOfTimesH0SIsAccepted"]], 25)
  checkEquals(result[["numberOfTimesH0SIsRejected"]], 0)
  checkEquals(result[["numberOfTimesFutilityDecidedAtLastLook"]], 0)
  checkEquals(result[["numberOfTimesTrialEndedAtLook"]],
          {out <- c(1, 9, 11, 4, 0); names(out) <- trialParameters$interimLookTime; out})
}


###################################################
### code chunk number 43: Case-F-unit-test
###################################################
test.caseF <- function() {
trialParameters <- list(minimumNumberOfEvents = 20,
                        minimumIncreaseInV = 0.2,
                        numberRecruitedEachYear = c(80, 120, 160, 160),
                        followupTime = 3,
                        adminCensoringTime = 7,
                        interimLookTime = c(1, 2, 3, 5, 7),
                        type1ErrorForResponse = 0.05,
                        type2ErrorForResponse = 0.01,
                        glrBoundarySidedness = "one", # one sided or two-sided
                        type1Error = 0.05,
                        type2Error = 0.10,
                        epsType1 = 1/3,
                        epsType2 = 1/3)
trueParameters <- list(p0 = 0.6,
                       p1 = 0.6,
                       pdiffHyp=0.3,
                       theta = list(
                         alpha = log(1.0),
                         beta = log(1.0),
                         gamma = log(1.0)),
                       baselineLambda = 0.35,
                       etaHyp = 0.25)
rngSeed <- 623312
sp23Design <- generateSP23Design(trueParameters, trialParameters)
print(sp23Design)
trialHistory <- exploreSP23Design(sp23Design, numberOfSimulations=25, rngSeed=rngSeed)
result <- analyzeSP23Design(sp23Design, trialHistory)$designSummary
  checkEquals(result[["numberOfTimesH0RIsRejectedAtFirstLook"]], 0)
  checkEquals(result[["numberOfTimesH0RIsRejected"]], 0)
  checkEquals(result[["numberOfTimesStoppedForFutility"]], 25)
  checkEquals(result[["numberOfTimesH0SIsAccepted"]], 0)
  checkEquals(result[["numberOfTimesH0SIsRejected"]], 0)
  checkEquals(result[["numberOfTimesFutilityDecidedAtLastLook"]], 0)
  checkEquals(result[["numberOfTimesTrialEndedAtLook"]],
          {out <- c(13, 8, 4, 0, 0); names(out) <- trialParameters$interimLookTime; out})
}


###################################################
### code chunk number 44: Case-G-unit-test
###################################################
test.caseG <- function() {
trialParameters <- list(minimumNumberOfEvents = 20,
                        minimumIncreaseInV = 0.2,
                        numberRecruitedEachYear = c(80, 120, 160, 160),
                        followupTime = 3,
                        adminCensoringTime = 7,
                        interimLookTime = c(1, 2, 3, 5, 7),
                        type1ErrorForResponse = 0.05,
                        type2ErrorForResponse = 0.01,
                        glrBoundarySidedness = "one", # one sided or two-sided
                        type1Error = 0.05,
                        type2Error = 0.10,
                        epsType1 = 1/3,
                        epsType2 = 1/3)
trueParameters <- list(p0 = 0.6,
                       p1 = 0.6,
                       pdiffHyp=0.3,
                       theta = list(
                         alpha = log(0.5),
                         beta = log(1.0),
                         gamma = log(1.0)),
                       baselineLambda = 0.35,
                       etaHyp = 0.25)
rngSeed <- 6232
sp23Design <- generateSP23Design(trueParameters, trialParameters)
print(sp23Design)
trialHistory <- exploreSP23Design(sp23Design, numberOfSimulations=25, rngSeed=rngSeed)
result <- analyzeSP23Design(sp23Design, trialHistory)$designSummary
  checkEquals(result[["numberOfTimesH0RIsRejectedAtFirstLook"]], 0)
  checkEquals(result[["numberOfTimesH0RIsRejected"]], 1)
  checkEquals(result[["numberOfTimesStoppedForFutility"]], 25)
  checkEquals(result[["numberOfTimesH0SIsAccepted"]], 1)
  checkEquals(result[["numberOfTimesH0SIsRejected"]], 0)
  checkEquals(result[["numberOfTimesFutilityDecidedAtLastLook"]], 0)
  checkEquals(result[["numberOfTimesTrialEndedAtLook"]],
          {out <- c(15, 8, 1, 1, 0); names(out) <- trialParameters$interimLookTime; out})
}


###################################################
### code chunk number 45: Run-1-unit-test
###################################################
test.run1 <- function() {
trialParameters <- list(minimumNumberOfEvents = 20,
                        minimumIncreaseInV = 0.2,
                        numberRecruitedEachYear = c(80, 120, 160, 160),
                        followupTime = 3,
                        adminCensoringTime = 7,
                        interimLookTime = c(1, 2, 3, 5, 7),
                        type1ErrorForResponse = 0.05,
                        type2ErrorForResponse = 0.01,
                        glrBoundarySidedness = "one", # one sided or two-sided
                        type1Error = 0.05,
                        type2Error = 0.10,
                        epsType1 = 1/3,
                        epsType2 = 1/3)
trueParameters <- list(p0 = 0.3,
                       p1 = 0.6,
                       pdiffHyp=0.3,
                       theta = list(
                         alpha = log(1.0),
                         beta = log(.75),
                         gamma = log(1.0)),
                       baselineLambda = 0.35,
                       etaHyp = 0.25)
rngSeed <- 5554231
sp23Design <- generateSP23Design(trueParameters, trialParameters)
print(sp23Design)
trialHistory <- exploreSP23Design(sp23Design, numberOfSimulations=25, rngSeed=rngSeed)
result <- analyzeSP23Design(sp23Design, trialHistory)$designSummary
  checkEquals(result[["numberOfTimesH0RIsRejectedAtFirstLook"]], 13)
  checkEquals(result[["numberOfTimesH0RIsRejected"]], 25)
  checkEquals(result[["numberOfTimesStoppedForFutility"]], 5)
  checkEquals(result[["numberOfTimesH0SIsAccepted"]], 5)
  checkEquals(result[["numberOfTimesH0SIsRejected"]], 20)
  checkEquals(result[["numberOfTimesFutilityDecidedAtLastLook"]], 4)
  checkEquals(result[["numberOfTimesTrialEndedAtLook"]],
          {out <- c(0, 3, 1, 5, 16); names(out) <- trialParameters$interimLookTime; out})
}


###################################################
### code chunk number 46: Run-2-unit-test
###################################################
test.run2 <- function() {
trialParameters <- list(minimumNumberOfEvents = 20,
                        minimumIncreaseInV = 0.2,
                        numberRecruitedEachYear = c(80, 120, 160, 160),
                        followupTime = 3,
                        adminCensoringTime = 7,
                        interimLookTime = c(1, 2, 3, 5, 7),
                        type1ErrorForResponse = 0.05,
                        type2ErrorForResponse = 0.01,
                        glrBoundarySidedness = "one", # one sided or two-sided
                        type1Error = 0.05,
                        type2Error = 0.10,
                        epsType1 = 1/3,
                        epsType2 = 1/3)
trueParameters <- list(p0 = 0.3,
                       p1 = 0.6,
                       pdiffHyp=0.3,
                       theta = list(
                         alpha = log(.75),
                         beta = log(.79),
                         gamma = log(1.0)),
                       baselineLambda = 0.35,
                       etaHyp = 0.25)
rngSeed <- 125352
sp23Design <- generateSP23Design(trueParameters, trialParameters)
print(sp23Design)
trialHistory <- exploreSP23Design(sp23Design, numberOfSimulations=25, rngSeed=rngSeed)
result <- analyzeSP23Design(sp23Design, trialHistory)$designSummary
  checkEquals(result[["numberOfTimesH0RIsRejectedAtFirstLook"]], 9)
  checkEquals(result[["numberOfTimesH0RIsRejected"]], 25)
  checkEquals(result[["numberOfTimesStoppedForFutility"]], 4)
  checkEquals(result[["numberOfTimesH0SIsAccepted"]], 4)
  checkEquals(result[["numberOfTimesH0SIsRejected"]], 21)
  checkEquals(result[["numberOfTimesFutilityDecidedAtLastLook"]], 2)
  checkEquals(result[["numberOfTimesTrialEndedAtLook"]],
          {out <- c(0, 4, 2, 10, 9); names(out) <- trialParameters$interimLookTime; out})
}


###################################################
### code chunk number 47: Run-3-unit-test
###################################################
test.run3 <- function() {
trialParameters <- list(minimumNumberOfEvents = 20,
                        minimumIncreaseInV = 0.2,
                        numberRecruitedEachYear = c(80, 120, 160, 160),
                        followupTime = 3,
                        adminCensoringTime = 7,
                        interimLookTime = c(1, 2, 3, 5, 7),
                        type1ErrorForResponse = 0.05,
                        type2ErrorForResponse = 0.01,
                        glrBoundarySidedness = "one", # one sided or two-sided
                        type1Error = 0.05,
                        type2Error = 0.10,
                        epsType1 = 1/3,
                        epsType2 = 1/3)
trueParameters <- list(p0 = 0.3,
                       p1 = 0.6,
                       pdiffHyp=0.3,
                       theta = list(
                         alpha = log(.5),
                         beta = log(.86),
                         gamma = log(1.0)),
                       baselineLambda = 0.35,
                       etaHyp = 0.25)
rngSeed <- 8669312
sp23Design <- generateSP23Design(trueParameters, trialParameters)
print(sp23Design)
trialHistory <- exploreSP23Design(sp23Design, numberOfSimulations=25, rngSeed=rngSeed)
result <- analyzeSP23Design(sp23Design, trialHistory)$designSummary
  checkEquals(result[["numberOfTimesH0RIsRejectedAtFirstLook"]], 11)
  checkEquals(result[["numberOfTimesH0RIsRejected"]], 25)
  checkEquals(result[["numberOfTimesStoppedForFutility"]], 1)
  checkEquals(result[["numberOfTimesH0SIsAccepted"]], 1)
  checkEquals(result[["numberOfTimesH0SIsRejected"]], 24)
  checkEquals(result[["numberOfTimesFutilityDecidedAtLastLook"]], 1)
  checkEquals(result[["numberOfTimesTrialEndedAtLook"]],
          {out <- c(0, 3, 3, 9, 10); names(out) <- trialParameters$interimLookTime; out})
}


###################################################
### code chunk number 48: Run-4-unit-test
###################################################
test.run4 <- function() {
trialParameters <- list(minimumNumberOfEvents = 20,
                        minimumIncreaseInV = 0.2,
                        numberRecruitedEachYear = c(80, 120, 160, 160),
                        followupTime = 3,
                        adminCensoringTime = 7,
                        interimLookTime = c(1, 2, 3, 5, 7),
                        type1ErrorForResponse = 0.05,
                        type2ErrorForResponse = 0.01,
                        glrBoundarySidedness = "one", # one sided or two-sided
                        type1Error = 0.05,
                        type2Error = 0.10,
                        epsType1 = 1/3,
                        epsType2 = 1/3)
trueParameters <- list(p0 = 0.3,
                       p1 = 0.6,
                       pdiffHyp=0.3,
                       theta = list(
                         alpha = log(.17),
                         beta = log(1.0),
                         gamma = log(1.0)),
                       baselineLambda = 0.35,
                       etaHyp = 0.25)
rngSeed <- 7325543
sp23Design <- generateSP23Design(trueParameters, trialParameters)
print(sp23Design)
trialHistory <- exploreSP23Design(sp23Design, numberOfSimulations=25, rngSeed=rngSeed)
result <- analyzeSP23Design(sp23Design, trialHistory)$designSummary
  checkEquals(result[["numberOfTimesH0RIsRejectedAtFirstLook"]], 16)
  checkEquals(result[["numberOfTimesH0RIsRejected"]], 25)
  checkEquals(result[["numberOfTimesStoppedForFutility"]], 6)
  checkEquals(result[["numberOfTimesH0SIsAccepted"]], 6)
  checkEquals(result[["numberOfTimesH0SIsRejected"]], 19)
  checkEquals(result[["numberOfTimesFutilityDecidedAtLastLook"]], 2)
  checkEquals(result[["numberOfTimesTrialEndedAtLook"]],
          {out <- c(0, 6, 2, 9, 8); names(out) <- trialParameters$interimLookTime; out})
}


###################################################
### code chunk number 49: Run-5-unit-test
###################################################
test.run5 <- function() {
trialParameters <- list(minimumNumberOfEvents = 20,
                        minimumIncreaseInV = 0.2,
                        numberRecruitedEachYear = c(80, 120, 160, 160),
                        followupTime = 3,
                        adminCensoringTime = 7,
                        interimLookTime = c(1, 2, 3, 5, 7),
                        type1ErrorForResponse = 0.05,
                        type2ErrorForResponse = 0.01,
                        glrBoundarySidedness = "one", # one sided or two-sided
                        type1Error = 0.05,
                        type2Error = 0.10,
                        epsType1 = 1/3,
                        epsType2 = 1/3)
trueParameters <- list(p0 = 0.3,
                       p1 = 0.6,
                       pdiffHyp=0.3,
                       theta = list(
                         alpha = log(.75),
                         beta = log(1.0),
                         gamma = log(0.61)),
                       baselineLambda = 0.35,
                       etaHyp = 0.25)
rngSeed <- 466553
sp23Design <- generateSP23Design(trueParameters, trialParameters)
print(sp23Design)
trialHistory <- exploreSP23Design(sp23Design, numberOfSimulations=25, rngSeed=rngSeed)
result <- analyzeSP23Design(sp23Design, trialHistory)$designSummary
  checkEquals(result[["numberOfTimesH0RIsRejectedAtFirstLook"]], 18)
  checkEquals(result[["numberOfTimesH0RIsRejected"]], 25)
  checkEquals(result[["numberOfTimesStoppedForFutility"]], 3)
  checkEquals(result[["numberOfTimesH0SIsAccepted"]], 3)
  checkEquals(result[["numberOfTimesH0SIsRejected"]], 22)
  checkEquals(result[["numberOfTimesFutilityDecidedAtLastLook"]], 2)
  checkEquals(result[["numberOfTimesTrialEndedAtLook"]],
          {out <- c(0, 4, 4, 6, 11); names(out) <- trialParameters$interimLookTime; out})
}


###################################################
### code chunk number 50: Run-6-unit-test
###################################################
test.run6 <- function() {
trialParameters <- list(minimumNumberOfEvents = 20,
                        minimumIncreaseInV = 0.2,
                        numberRecruitedEachYear = c(80, 120, 160, 160),
                        followupTime = 3,
                        adminCensoringTime = 7,
                        interimLookTime = c(1, 2, 3, 5, 7),
                        type1ErrorForResponse = 0.05,
                        type2ErrorForResponse = 0.01,
                        glrBoundarySidedness = "one", # one sided or two-sided
                        type1Error = 0.05,
                        type2Error = 0.10,
                        epsType1 = 1/3,
                        epsType2 = 1/3)
trueParameters <- list(p0 = 0.3,
                       p1 = 0.6,
                       pdiffHyp=0.3,
                       theta = list(
                         alpha = log(.5),
                         beta = log(1.0),
                         gamma = log(0.67)),
                       baselineLambda = 0.35,
                       etaHyp = 0.25)
rngSeed <- 1875543
sp23Design <- generateSP23Design(trueParameters, trialParameters)
print(sp23Design)
trialHistory <- exploreSP23Design(sp23Design, numberOfSimulations=25, rngSeed=rngSeed)
result <- analyzeSP23Design(sp23Design, trialHistory)$designSummary
  checkEquals(result[["numberOfTimesH0RIsRejectedAtFirstLook"]], 17)
  checkEquals(result[["numberOfTimesH0RIsRejected"]], 25)
  checkEquals(result[["numberOfTimesStoppedForFutility"]], 2)
  checkEquals(result[["numberOfTimesH0SIsAccepted"]], 2)
  checkEquals(result[["numberOfTimesH0SIsRejected"]], 23)
  checkEquals(result[["numberOfTimesFutilityDecidedAtLastLook"]], 1)
  checkEquals(result[["numberOfTimesTrialEndedAtLook"]],
          {out <- c(0, 4, 5, 9, 7); names(out) <- trialParameters$interimLookTime; out})
}


###################################################
### code chunk number 51: Run-7-unit-test
###################################################
test.run7 <- function() {
trialParameters <- list(minimumNumberOfEvents = 20,
                        minimumIncreaseInV = 0.2,
                        numberRecruitedEachYear = c(80, 120, 160, 160),
                        followupTime = 3,
                        adminCensoringTime = 7,
                        interimLookTime = c(1, 2, 3, 5, 7),
                        type1ErrorForResponse = 0.05,
                        type2ErrorForResponse = 0.01,
                        glrBoundarySidedness = "one", # one sided or two-sided
                        type1Error = 0.05,
                        type2Error = 0.10,
                        epsType1 = 1/3,
                        epsType2 = 1/3)
trueParameters <- list(p0 = 0.3,
                       p1 = 0.5,
                       pdiffHyp=0.3,
                       theta = list(
                         alpha = log(.75),
                         beta = log(1.0),
                         gamma = log(0.47)),
                       baselineLambda = 0.35,
                       etaHyp = 0.25)
rngSeed <- 1093
sp23Design <- generateSP23Design(trueParameters, trialParameters)
print(sp23Design)
trialHistory <- exploreSP23Design(sp23Design, numberOfSimulations=25, rngSeed=rngSeed)
result <- analyzeSP23Design(sp23Design, trialHistory)$designSummary
  checkEquals(result[["numberOfTimesH0RIsRejectedAtFirstLook"]], 5)
  checkEquals(result[["numberOfTimesH0RIsRejected"]], 25)
  checkEquals(result[["numberOfTimesStoppedForFutility"]], 3)
  checkEquals(result[["numberOfTimesH0SIsAccepted"]], 3)
  checkEquals(result[["numberOfTimesH0SIsRejected"]], 22)
  checkEquals(result[["numberOfTimesFutilityDecidedAtLastLook"]], 1)
  checkEquals(result[["numberOfTimesTrialEndedAtLook"]],
          {out <- c(0, 3, 3, 5, 14); names(out) <- trialParameters$interimLookTime; out})
}


###################################################
### code chunk number 52: Run-8-unit-test
###################################################
test.run8 <- function() {
trialParameters <- list(minimumNumberOfEvents = 20,
                        minimumIncreaseInV = 0.2,
                        numberRecruitedEachYear = c(80, 120, 160, 160),
                        followupTime = 3,
                        adminCensoringTime = 7,
                        interimLookTime = c(1, 2, 3, 5, 7),
                        type1ErrorForResponse = 0.05,
                        type2ErrorForResponse = 0.01,
                        glrBoundarySidedness = "one", # one sided or two-sided
                        type1Error = 0.05,
                        type2Error = 0.10,
                        epsType1 = 1/3,
                        epsType2 = 1/3)
trueParameters <- list(p0 = 0.3,
                       p1 = 0.5,
                       pdiffHyp=0.3,
                       theta = list(
                         alpha = 0,
                         beta = log(.75),
                         gamma = 0),
                       baselineLambda = 0.35,
                       etaHyp = 0.25)
rngSeed <- 117
sp23Design <- generateSP23Design(trueParameters, trialParameters)
print(sp23Design)
trialHistory <- exploreSP23Design(sp23Design, numberOfSimulations=25, rngSeed=rngSeed)
result <- analyzeSP23Design(sp23Design, trialHistory)$designSummary
  checkEquals(result[["numberOfTimesH0RIsRejectedAtFirstLook"]], 9)
  checkEquals(result[["numberOfTimesH0RIsRejected"]], 23)
  checkEquals(result[["numberOfTimesStoppedForFutility"]], 8)
  checkEquals(result[["numberOfTimesH0SIsAccepted"]], 6)
  checkEquals(result[["numberOfTimesH0SIsRejected"]], 17)
  checkEquals(result[["numberOfTimesFutilityDecidedAtLastLook"]], 3)
  checkEquals(result[["numberOfTimesTrialEndedAtLook"]],
          {out <- c(0, 5, 7, 7, 6); names(out) <- trialParameters$interimLookTime; out})
}


###################################################
### code chunk number 53: Run-9-unit-test
###################################################
test.run9 <- function() {
trialParameters <- list(minimumNumberOfEvents = 20,
                        minimumIncreaseInV = 0.2,
                        numberRecruitedEachYear = c(80, 120, 160, 160),
                        followupTime = 3,
                        adminCensoringTime = 7,
                        interimLookTime = c(1, 2, 3, 5, 7),
                        type1ErrorForResponse = 0.05,
                        type2ErrorForResponse = 0.01,
                        glrBoundarySidedness = "one", # one sided or two-sided
                        type1Error = 0.05,
                        type2Error = 0.10,
                        epsType1 = 1/3,
                        epsType2 = 1/3)
trueParameters <- list(p0 = 0.2,
                       p1 = 0.5,
                       pdiffHyp=0.3,
                       theta = list(
                         alpha = 0,
                         beta = log(.75),
                         gamma = 0),
                       baselineLambda = 0.35,
                       etaHyp = 0.25)
rngSeed <- 9872831
sp23Design <- generateSP23Design(trueParameters, trialParameters)
print(sp23Design)
trialHistory <- exploreSP23Design(sp23Design, numberOfSimulations=25, rngSeed=rngSeed)
result <- analyzeSP23Design(sp23Design, trialHistory)$designSummary
  checkEquals(result[["numberOfTimesH0RIsRejectedAtFirstLook"]], 16)
  checkEquals(result[["numberOfTimesH0RIsRejected"]], 25)
  checkEquals(result[["numberOfTimesStoppedForFutility"]], 3)
  checkEquals(result[["numberOfTimesH0SIsAccepted"]], 3)
  checkEquals(result[["numberOfTimesH0SIsRejected"]], 22)
  checkEquals(result[["numberOfTimesFutilityDecidedAtLastLook"]], 3)
  checkEquals(result[["numberOfTimesTrialEndedAtLook"]],
          {out <- c(0, 1, 4, 9, 11); names(out) <- trialParameters$interimLookTime; out})
}


