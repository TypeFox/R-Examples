library(sp23design)
trueParameters <- list(p0 = 0.3,
                       p1 = 0.6,
                       pdiffHyp=0.3,
                       theta = list(
                         alpha = log(0.2),
                         beta = log(1.56),
                         gamma = log(1.0)),
                       baselineLambda = 0.35,
                       etaHyp = 0.25)

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

sp23Design <- generateSP23Design(trueParameters, trialParameters)
trialHistory <- exploreSP23Design(sp23Design, numberOfSimulations=1000, rngSeed=283119)
result <- analyzeSP23Design(sp23Design, trialHistory)$designSummary

cat("numberOfTimesH0RIsRejectedAtFirstLook", result[["numberOfTimesH0RIsRejectedAtFirstLook"]], "\n")
cat("numberOfTimesH0RIsRejected", result[["numberOfTimesH0RIsRejected"]], "\n")
cat("numberOfTimesStoppedForFutility", result[["numberOfTimesStoppedForFutility"]], "\n")
cat("numberOfTimesH0SIsAccepted", result[["numberOfTimesH0SIsAccepted"]], "\n")
cat("numberOfTimesH0SIsRejected", result[["numberOfTimesH0SIsRejected"]], "\n")
cat("numberOfTimesFutilityDecidedAtLastLook", result[["numberOfTimesFutilityDecidedAtLastLook"]], "\n")
cat("numberOfTimesTrialEndedAtLook", result[["numberOfTimesTrialEndedAtLook"]], "\n")            
cat("avgExitTime", result[["avgExitTime"]], "\n")
     

