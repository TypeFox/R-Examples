pfnorm <- function(data){
  
  s <- data.frame(ID=data[,1],
                  Age=data[,2],
                  Gender=data[,3],
                  NumberofSteps=data[,4],
                  Numberoffullstands=data[,5],
                  UpGo8ftmean=data[,6],
                  armcurlR_Lmean=data[,7])
  
  s$agecat <- rep(NA,nrow(s))
  s$agecat[s$Age >= 60 & s$Age < 65] <- 1
  s$agecat[s$Age >= 65 & s$Age < 70] <- 2
  s$agecat[s$Age >= 70 & s$Age < 75] <- 3
  s$agecat[s$Age >= 75 & s$Age < 80] <- 4
  s$agecat[s$Age >= 80 & s$Age < 85] <- 5
  s$agecat[s$Age >= 85 & s$Age < 90] <- 6
  s$agecat[s$Age >= 90 & s$Age < 94] <- 7
  
  s$n.stand[s$Gender == 'F' & s$agecat == 1 & s$Numberoffullstands >= 12 & s$Numberoffullstands <= 17] <- 'Meets the norm'
  s$n.stand[s$Gender == 'F' & s$agecat == 1 & s$Numberoffullstands < 12] <- 'Lower than the norm'
  s$n.stand[s$Gender == 'F' & s$agecat == 1 & s$Numberoffullstands > 17] <- 'Higher than the norm'
  
  s$n.stand[s$Gender == 'F' & s$agecat == 2 & s$Numberoffullstands >= 11 & s$Numberoffullstands <= 16] <- 'Meets the norm'
  s$n.stand[s$Gender == 'F' & s$agecat == 2 & s$Numberoffullstands < 11] <- 'Lower than the norm'
  s$n.stand[s$Gender == 'F' & s$agecat == 2 & s$Numberoffullstands > 16] <- 'Higher than the norm'
  
  s$n.stand[s$Gender == 'F' & s$agecat == 3 & s$Numberoffullstands >= 10 & s$Numberoffullstands <= 15] <- 'Meets the norm'
  s$n.stand[s$Gender == 'F' & s$agecat == 3 & s$Numberoffullstands < 10] <- 'Lower than the norm'
  s$n.stand[s$Gender == 'F' & s$agecat == 3 & s$Numberoffullstands > 15] <- 'Higher than the norm'
  
  s$n.stand[s$Gender == 'F' & s$agecat == 4 & s$Numberoffullstands >=10 & s$Numberoffullstands <= 15] <- 'Meets the norm'
  s$n.stand[s$Gender == 'F' & s$agecat == 4 & s$Numberoffullstands <10] <- 'Lower than the norm'
  s$n.stand[s$Gender == 'F' & s$agecat == 4 & s$Numberoffullstands > 15] <- 'Higher than the norm'
  
  s$n.stand[s$Gender == 'F' & s$agecat == 5 & s$Numberoffullstands >=9 & s$Numberoffullstands <= 14] <- 'Meets the norm'
  s$n.stand[s$Gender == 'F' & s$agecat == 5 & s$Numberoffullstands <9] <- 'Lower than the norm'
  s$n.stand[s$Gender == 'F' & s$agecat == 5 & s$Numberoffullstands > 14] <- 'Higher than the norm'
  
  s$n.stand[s$Gender == 'F' & s$agecat == 6 & s$Numberoffullstands >=8 & s$Numberoffullstands <= 13] <- 'Meets the norm'
  s$n.stand[s$Gender == 'F' & s$agecat == 6 & s$Numberoffullstands <8] <- 'Lower than the norm'
  s$n.stand[s$Gender == 'F' & s$agecat == 6 & s$Numberoffullstands > 13] <- 'Higher than the norm'
  
  s$n.stand[s$Gender == 'F' & s$agecat == 7 & s$Numberoffullstands >=4 & s$Numberoffullstands <= 11] <- 'Meets the norm'
  s$n.stand[s$Gender == 'F' & s$agecat == 7 & s$Numberoffullstands <4] <- 'Lower than the norm'
  s$n.stand[s$Gender == 'F' & s$agecat == 7 & s$Numberoffullstands > 11] <- 'Higher than the norm'
  
  s$n.stand[s$Gender == 'M' & s$agecat == 1 & s$Numberoffullstands >= 14 & s$Numberoffullstands <= 19] <- 'Meets the norm'
  s$n.stand[s$Gender == 'M' & s$agecat == 1 & s$Numberoffullstands < 14] <- 'Lower than the norm'
  s$n.stand[s$Gender == 'M' & s$agecat == 1 & s$Numberoffullstands > 19] <- 'Higher than the norm'
  
  s$n.stand[s$Gender == 'M' & s$agecat == 2 & s$Numberoffullstands >= 12 & s$Numberoffullstands <= 18] <- 'Meets the norm'
  s$n.stand[s$Gender == 'M' & s$agecat == 2 & s$Numberoffullstands < 12] <- 'Lower than the norm'
  s$n.stand[s$Gender == 'M' & s$agecat == 2 & s$Numberoffullstands > 18] <- 'Higher than the norm'
  
  s$n.stand[s$Gender == 'M' & s$agecat == 3 & s$Numberoffullstands >= 12 & s$Numberoffullstands <= 17] <- 'Meets the norm'
  s$n.stand[s$Gender == 'M' & s$agecat == 3 & s$Numberoffullstands < 12] <- 'Lower than the norm'
  s$n.stand[s$Gender == 'M' & s$agecat == 3 & s$Numberoffullstands > 17] <- 'Higher than the norm'
  
  s$n.stand[s$Gender == 'M' & s$agecat == 4 & s$Numberoffullstands >=11 & s$Numberoffullstands <= 17] <- 'Meets the norm'
  s$n.stand[s$Gender == 'M' & s$agecat == 4 & s$Numberoffullstands <11] <- 'Lower than the norm'
  s$n.stand[s$Gender == 'M' & s$agecat == 4 & s$Numberoffullstands > 17] <- 'Higher than the norm'
  
  s$n.stand[s$Gender == 'M' & s$agecat == 5 & s$Numberoffullstands >=10 & s$Numberoffullstands <= 15] <- 'Meets the norm'
  s$n.stand[s$Gender == 'M' & s$agecat == 5 & s$Numberoffullstands <10] <- 'Lower than the norm'
  s$n.stand[s$Gender == 'M' & s$agecat == 5 & s$Numberoffullstands > 15] <- 'Higher than the norm'
  
  s$n.stand[s$Gender == 'M' & s$agecat == 6 & s$Numberoffullstands >=8 & s$Numberoffullstands <= 14] <- 'Meets the norm'
  s$n.stand[s$Gender == 'M' & s$agecat == 6 & s$Numberoffullstands <8] <- 'Lower than the norm'
  s$n.stand[s$Gender == 'M' & s$agecat == 6 & s$Numberoffullstands > 14] <- 'Higher than the norm'
  
  s$n.stand[s$Gender == 'M' & s$agecat == 7 & s$Numberoffullstands >=7 & s$Numberoffullstands <= 12] <- 'Meets the norm'
  s$n.stand[s$Gender == 'M' & s$agecat == 7 & s$Numberoffullstands <7] <- 'Lower than the norm'
  s$n.stand[s$Gender == 'M' & s$agecat == 7 & s$Numberoffullstands > 12] <- 'Higher than the norm'
  
  s$n.armcurl[s$Gender == 'F' & s$agecat == 1 & s$armcurlR_Lmean >= 13 & s$armcurlR_Lmean <= 19] <- 'Meets the norm'
  s$n.armcurl[s$Gender == 'F' & s$agecat == 1 & s$armcurlR_Lmean < 13] <- 'Lower than the norm'
  s$n.armcurl[s$Gender == 'F' & s$agecat == 1 & s$armcurlR_Lmean > 19] <- 'Higher than the norm'
  
  s$n.armcurl[s$Gender == 'F' & s$agecat == 2 & s$armcurlR_Lmean >= 12 & s$armcurlR_Lmean <= 18] <- 'Meets the norm'
  s$n.armcurl[s$Gender == 'F' & s$agecat == 2 & s$armcurlR_Lmean < 12] <- 'Lower than the norm'
  s$n.armcurl[s$Gender == 'F' & s$agecat == 2 & s$armcurlR_Lmean > 18] <- 'Higher than the norm'
  
  s$n.armcurl[s$Gender == 'F' & s$agecat == 3 & s$armcurlR_Lmean >= 12 & s$armcurlR_Lmean <= 17] <- 'Meets the norm'
  s$n.armcurl[s$Gender == 'F' & s$agecat == 3 & s$armcurlR_Lmean < 12] <- 'Lower than the norm'
  s$n.armcurl[s$Gender == 'F' & s$agecat == 3 & s$armcurlR_Lmean > 17] <- 'Higher than the norm'
  
  s$n.armcurl[s$Gender == 'F' & s$agecat == 4 & s$armcurlR_Lmean >=11 & s$armcurlR_Lmean <= 17] <- 'Meets the norm'
  s$n.armcurl[s$Gender == 'F' & s$agecat == 4 & s$armcurlR_Lmean <11] <- 'Lower than the norm'
  s$n.armcurl[s$Gender == 'F' & s$agecat == 4 & s$armcurlR_Lmean > 17] <- 'Higher than the norm'
  
  s$n.armcurl[s$Gender == 'F' & s$agecat == 5 & s$armcurlR_Lmean >=10 & s$armcurlR_Lmean <= 16] <- 'Meets the norm'
  s$n.armcurl[s$Gender == 'F' & s$agecat == 5 & s$armcurlR_Lmean <10] <- 'Lower than the norm'
  s$n.armcurl[s$Gender == 'F' & s$agecat == 5 & s$armcurlR_Lmean > 16] <- 'Higher than the norm'
  
  s$n.armcurl[s$Gender == 'F' & s$agecat == 6 & s$armcurlR_Lmean >=10 & s$armcurlR_Lmean <= 15] <- 'Meets the norm'
  s$n.armcurl[s$Gender == 'F' & s$agecat == 6 & s$armcurlR_Lmean <10] <- 'Lower than the norm'
  s$n.armcurl[s$Gender == 'F' & s$agecat == 6 & s$armcurlR_Lmean > 15] <- 'Higher than the norm'
  
  s$n.armcurl[s$Gender == 'F' & s$agecat == 7 & s$armcurlR_Lmean >=8 & s$armcurlR_Lmean <= 13] <- 'Meets the norm'
  s$n.armcurl[s$Gender == 'F' & s$agecat == 7 & s$armcurlR_Lmean <8] <- 'Lower than the norm'
  s$n.armcurl[s$Gender == 'F' & s$agecat == 7 & s$armcurlR_Lmean > 13] <- 'Higher than the norm'
  
  s$n.armcurl[s$Gender == 'M' & s$agecat == 1 & s$armcurlR_Lmean >= 16 & s$armcurlR_Lmean <= 22] <- 'Meets the norm'
  s$n.armcurl[s$Gender == 'M' & s$agecat == 1 & s$armcurlR_Lmean < 16] <- 'Lower than the norm'
  s$n.armcurl[s$Gender == 'M' & s$agecat == 1 & s$armcurlR_Lmean > 22] <- 'Higher than the norm'
  
  s$n.armcurl[s$Gender == 'M' & s$agecat == 2 & s$armcurlR_Lmean >= 15 & s$armcurlR_Lmean <= 21] <- 'Meets the norm'
  s$n.armcurl[s$Gender == 'M' & s$agecat == 2 & s$armcurlR_Lmean < 15] <- 'Lower than the norm'
  s$n.armcurl[s$Gender == 'M' & s$agecat == 2 & s$armcurlR_Lmean > 21] <- 'Higher than the norm'
  
  s$n.armcurl[s$Gender == 'M' & s$agecat == 3 & s$armcurlR_Lmean >= 14 & s$armcurlR_Lmean <= 21] <- 'Meets the norm'
  s$n.armcurl[s$Gender == 'M' & s$agecat == 3 & s$armcurlR_Lmean < 14] <- 'Lower than the norm'
  s$n.armcurl[s$Gender == 'M' & s$agecat == 3 & s$armcurlR_Lmean > 21] <- 'Higher than the norm'
  
  s$n.armcurl[s$Gender == 'M' & s$agecat == 4 & s$armcurlR_Lmean >=13 & s$armcurlR_Lmean <= 19] <- 'Meets the norm'
  s$n.armcurl[s$Gender == 'M' & s$agecat == 4 & s$armcurlR_Lmean <13] <- 'Lower than the norm'
  s$n.armcurl[s$Gender == 'M' & s$agecat == 4 & s$armcurlR_Lmean > 19] <- 'Higher than the norm'
  
  s$n.armcurl[s$Gender == 'M' & s$agecat == 5 & s$armcurlR_Lmean >=13 & s$armcurlR_Lmean <= 19] <- 'Meets the norm'
  s$n.armcurl[s$Gender == 'M' & s$agecat == 5 & s$armcurlR_Lmean <13] <- 'Lower than the norm'
  s$n.armcurl[s$Gender == 'M' & s$agecat == 5 & s$armcurlR_Lmean > 19] <- 'Higher than the norm'
  
  s$n.armcurl[s$Gender == 'M' & s$agecat == 6 & s$armcurlR_Lmean >=11 & s$armcurlR_Lmean <= 17] <- 'Meets the norm'
  s$n.armcurl[s$Gender == 'M' & s$agecat == 6 & s$armcurlR_Lmean <11] <- 'Lower than the norm'
  s$n.armcurl[s$Gender == 'M' & s$agecat == 6 & s$armcurlR_Lmean > 17] <- 'Higher than the norm'
  
  s$n.armcurl[s$Gender == 'M' & s$agecat == 7 & s$armcurlR_Lmean >=10 & s$armcurlR_Lmean <= 14] <- 'Meets the norm'
  s$n.armcurl[s$Gender == 'M' & s$agecat == 7 & s$armcurlR_Lmean <10] <- 'Lower than the norm'
  s$n.armcurl[s$Gender == 'M' & s$agecat == 7 & s$armcurlR_Lmean > 14] <- 'Higher than the norm'
  
  s$n.steps[s$Gender == 'F' & s$agecat == 1 & s$NumberofSteps >= 75 & s$NumberofSteps <= 107] <- 'Meets the norm'
  s$n.steps[s$Gender == 'F' & s$agecat == 1 & s$NumberofSteps < 75] <- 'Lower than the norm'
  s$n.steps[s$Gender == 'F' & s$agecat == 1 & s$NumberofSteps > 107] <- 'Higher than the norm'
  
  s$n.steps[s$Gender == 'F' & s$agecat == 2 & s$NumberofSteps >= 73 & s$NumberofSteps <= 107] <- 'Meets the norm'
  s$n.steps[s$Gender == 'F' & s$agecat == 2 & s$NumberofSteps < 73] <- 'Lower than the norm'
  s$n.steps[s$Gender == 'F' & s$agecat == 2 & s$NumberofSteps > 107] <- 'Higher than the norm'
  
  s$n.steps[s$Gender == 'F' & s$agecat == 3 & s$NumberofSteps >= 68 & s$NumberofSteps <= 101] <- 'Meets the norm'
  s$n.steps[s$Gender == 'F' & s$agecat == 3 & s$NumberofSteps < 68] <- 'Lower than the norm'
  s$n.steps[s$Gender == 'F' & s$agecat == 3 & s$NumberofSteps > 101] <- 'Higher than the norm'
  
  s$n.steps[s$Gender == 'F' & s$agecat == 4 & s$NumberofSteps >=68 & s$NumberofSteps <= 100] <- 'Meets the norm'
  s$n.steps[s$Gender == 'F' & s$agecat == 4 & s$NumberofSteps <68] <- 'Lower than the norm'
  s$n.steps[s$Gender == 'F' & s$agecat == 4 & s$NumberofSteps > 100] <- 'Higher than the norm'
  
  s$n.steps[s$Gender == 'F' & s$agecat == 5 & s$NumberofSteps >=60 & s$NumberofSteps <= 90] <- 'Meets the norm'
  s$n.steps[s$Gender == 'F' & s$agecat == 5 & s$NumberofSteps <60] <- 'Lower than the norm'
  s$n.steps[s$Gender == 'F' & s$agecat == 5 & s$NumberofSteps > 90] <- 'Higher than the norm'
  
  s$n.steps[s$Gender == 'F' & s$agecat == 6 & s$NumberofSteps >=55 & s$NumberofSteps <= 85] <- 'Meets the norm'
  s$n.steps[s$Gender == 'F' & s$agecat == 6 & s$NumberofSteps <55] <- 'Lower than the norm'
  s$n.steps[s$Gender == 'F' & s$agecat == 6 & s$NumberofSteps > 85] <- 'Higher than the norm'
  
  s$n.steps[s$Gender == 'F' & s$agecat == 7 & s$NumberofSteps >=44 & s$NumberofSteps <= 72] <- 'Meets the norm'
  s$n.steps[s$Gender == 'F' & s$agecat == 7 & s$NumberofSteps <44] <- 'Lower than the norm'
  s$n.steps[s$Gender == 'F' & s$agecat == 7 & s$NumberofSteps > 72] <- 'Higher than the norm'
  
  s$n.steps[s$Gender == 'M' & s$agecat == 1 & s$NumberofSteps >= 87 & s$NumberofSteps <= 115] <- 'Meets the norm'
  s$n.steps[s$Gender == 'M' & s$agecat == 1 & s$NumberofSteps < 87] <- 'Lower than the norm'
  s$n.steps[s$Gender == 'M' & s$agecat == 1 & s$NumberofSteps > 115] <- 'Higher than the norm'
  
  s$n.steps[s$Gender == 'M' & s$agecat == 2 & s$NumberofSteps >= 86 & s$NumberofSteps <= 116] <- 'Meets the norm'
  s$n.steps[s$Gender == 'M' & s$agecat == 2 & s$NumberofSteps < 86] <- 'Lower than the norm'
  s$n.steps[s$Gender == 'M' & s$agecat == 2 & s$NumberofSteps > 116] <- 'Higher than the norm'
  
  s$n.steps[s$Gender == 'M' & s$agecat == 3 & s$NumberofSteps >= 80 & s$NumberofSteps <= 110] <- 'Meets the norm'
  s$n.steps[s$Gender == 'M' & s$agecat == 3 & s$NumberofSteps < 80] <- 'Lower than the norm'
  s$n.steps[s$Gender == 'M' & s$agecat == 3 & s$NumberofSteps > 110] <- 'Higher than the norm'
  
  s$n.steps[s$Gender == 'M' & s$agecat == 4 & s$NumberofSteps >=73 & s$NumberofSteps <= 109] <- 'Meets the norm'
  s$n.steps[s$Gender == 'M' & s$agecat == 4 & s$NumberofSteps <73] <- 'Lower than the norm'
  s$n.steps[s$Gender == 'M' & s$agecat == 4 & s$NumberofSteps > 109] <- 'Higher than the norm'
  
  s$n.steps[s$Gender == 'M' & s$agecat == 5 & s$NumberofSteps >=71 & s$NumberofSteps <= 103] <- 'Meets the norm'
  s$n.steps[s$Gender == 'M' & s$agecat == 5 & s$NumberofSteps <71] <- 'Lower than the norm'
  s$n.steps[s$Gender == 'M' & s$agecat == 5 & s$NumberofSteps > 103] <- 'Higher than the norm'
  
  s$n.steps[s$Gender == 'M' & s$agecat == 6 & s$NumberofSteps >=59 & s$NumberofSteps <= 91] <- 'Meets the norm'
  s$n.steps[s$Gender == 'M' & s$agecat == 6 & s$NumberofSteps <59] <- 'Lower than the norm'
  s$n.steps[s$Gender == 'M' & s$agecat == 6 & s$NumberofSteps > 91] <- 'Higher than the norm'
  
  s$n.steps[s$Gender == 'M' & s$agecat == 7 & s$NumberofSteps >=52 & s$NumberofSteps <= 86] <- 'Meets the norm'
  s$n.steps[s$Gender == 'M' & s$agecat == 7 & s$NumberofSteps <52] <- 'Lower than the norm'
  s$n.steps[s$Gender == 'M' & s$agecat == 7 & s$NumberofSteps > 86] <- 'Higher than the norm'
  
  s$n.upgo[s$Gender == 'F' & s$agecat == 1 & s$UpGo8ftmean >= 4.4 & s$UpGo8ftmean <= 6] <- 'Meets the norm'
  s$n.upgo[s$Gender == 'F' & s$agecat == 1 & s$UpGo8ftmean < 4.4] <- 'Higher than the norm'
  s$n.upgo[s$Gender == 'F' & s$agecat == 1 & s$UpGo8ftmean > 6] <- 'Lower than the norm'
  
  s$n.upgo[s$Gender == 'F' & s$agecat == 2 & s$UpGo8ftmean >= 4.8 & s$UpGo8ftmean <= 6.4] <- 'Meets the norm'
  s$n.upgo[s$Gender == 'F' & s$agecat == 2 & s$UpGo8ftmean < 4.8] <- 'Higher than the norm'
  s$n.upgo[s$Gender == 'F' & s$agecat == 2 & s$UpGo8ftmean > 6.4] <- 'Lower than the norm'
  
  s$n.upgo[s$Gender == 'F' & s$agecat == 3 & s$UpGo8ftmean >= 4.9 & s$UpGo8ftmean <= 7.1] <- 'Meets the norm'
  s$n.upgo[s$Gender == 'F' & s$agecat == 3 & s$UpGo8ftmean < 4.9] <- 'Higher than the norm'
  s$n.upgo[s$Gender == 'F' & s$agecat == 3 & s$UpGo8ftmean > 7.1] <- 'Lower than the norm'
  
  s$n.upgo[s$Gender == 'F' & s$agecat == 4 & s$UpGo8ftmean >=5.2 & s$UpGo8ftmean <= 7.4] <- 'Meets the norm'
  s$n.upgo[s$Gender == 'F' & s$agecat == 4 & s$UpGo8ftmean <5.2] <- 'Higher than the norm'
  s$n.upgo[s$Gender == 'F' & s$agecat == 4 & s$UpGo8ftmean > 7.4] <- 'Lower than the norm'
  
  s$n.upgo[s$Gender == 'F' & s$agecat == 5 & s$UpGo8ftmean >=5.7 & s$UpGo8ftmean <= 8.7] <- 'Meets the norm'
  s$n.upgo[s$Gender == 'F' & s$agecat == 5 & s$UpGo8ftmean <5.7] <- 'Higher than the norm'
  s$n.upgo[s$Gender == 'F' & s$agecat == 5 & s$UpGo8ftmean > 8.7] <- 'Lower than the norm'
  
  s$n.upgo[s$Gender == 'F' & s$agecat == 6 & s$UpGo8ftmean >=6.2 & s$UpGo8ftmean <= 9.6] <- 'Meets the norm'
  s$n.upgo[s$Gender == 'F' & s$agecat == 6 & s$UpGo8ftmean <6.2] <- 'Higher than the norm'
  s$n.upgo[s$Gender == 'F' & s$agecat == 6 & s$UpGo8ftmean > 9.6] <- 'Lower than the norm'
  
  s$n.upgo[s$Gender == 'F' & s$agecat == 7 & s$UpGo8ftmean >=7.3 & s$UpGo8ftmean <= 11.5] <- 'Meets the norm'
  s$n.upgo[s$Gender == 'F' & s$agecat == 7 & s$UpGo8ftmean <7.3] <- 'Higher than the norm'
  s$n.upgo[s$Gender == 'F' & s$agecat == 7 & s$UpGo8ftmean > 11.5] <- 'Lower than the norm'
  
  s$n.upgo[s$Gender == 'M' & s$agecat == 1 & s$UpGo8ftmean >= 3.8 & s$UpGo8ftmean <= 5.6] <- 'Meets the norm'
  s$n.upgo[s$Gender == 'M' & s$agecat == 1 & s$UpGo8ftmean < 3.8] <- 'Higher than the norm'
  s$n.upgo[s$Gender == 'M' & s$agecat == 1 & s$UpGo8ftmean > 5.6] <- 'Lower than the norm'
  
  s$n.upgo[s$Gender == 'M' & s$agecat == 2 & s$UpGo8ftmean >= 4.3 & s$UpGo8ftmean <= 5.9] <- 'Meets the norm'
  s$n.upgo[s$Gender == 'M' & s$agecat == 2 & s$UpGo8ftmean < 4.3] <- 'Higher than the norm'
  s$n.upgo[s$Gender == 'M' & s$agecat == 2 & s$UpGo8ftmean > 5.9] <- 'Lower than the norm'
  
  s$n.upgo[s$Gender == 'M' & s$agecat == 3 & s$UpGo8ftmean >= 4.4 & s$UpGo8ftmean <= 6.2] <- 'Meets the norm'
  s$n.upgo[s$Gender == 'M' & s$agecat == 3 & s$UpGo8ftmean < 4.4] <- 'Higher than the norm'
  s$n.upgo[s$Gender == 'M' & s$agecat == 3 & s$UpGo8ftmean > 6.2] <- 'Lower than the norm'
  
  s$n.upgo[s$Gender == 'M' & s$agecat == 4 & s$UpGo8ftmean >=4.6 & s$UpGo8ftmean <= 7.2] <- 'Meets the norm'
  s$n.upgo[s$Gender == 'M' & s$agecat == 4 & s$UpGo8ftmean <4.6] <- 'Higher than the norm'
  s$n.upgo[s$Gender == 'M' & s$agecat == 4 & s$UpGo8ftmean > 7.2] <- 'Lower than the norm'
  
  s$n.upgo[s$Gender == 'M' & s$agecat == 5 & s$UpGo8ftmean >=5.2 & s$UpGo8ftmean <= 7.6] <- 'Meets the norm'
  s$n.upgo[s$Gender == 'M' & s$agecat == 5 & s$UpGo8ftmean <5.2] <- 'Higher than the norm'
  s$n.upgo[s$Gender == 'M' & s$agecat == 5 & s$UpGo8ftmean > 7.6] <- 'Lower than the norm'
  
  s$n.upgo[s$Gender == 'M' & s$agecat == 6 & s$UpGo8ftmean >=5.5 & s$UpGo8ftmean <= 8.9] <- 'Meets the norm'
  s$n.upgo[s$Gender == 'M' & s$agecat == 6 & s$UpGo8ftmean <5.5] <- 'Higher than the norm'
  s$n.upgo[s$Gender == 'M' & s$agecat == 6 & s$UpGo8ftmean > 8.9] <- 'Lower than the norm'
  
  s$n.upgo[s$Gender == 'M' & s$agecat == 7 & s$UpGo8ftmean >=6.2 & s$UpGo8ftmean <= 10] <- 'Meets the norm'
  s$n.upgo[s$Gender == 'M' & s$agecat == 7 & s$UpGo8ftmean <6.2] <- 'Higher than the norm'
  s$n.upgo[s$Gender == 'M' & s$agecat == 7 & s$UpGo8ftmean > 10] <- 'Lower than the norm'
  
  norm <- data.frame(ID = s[,1],
                     Age = s[,2],
                     Gender = s[,3],
                     Steps = s[,4],
                     Stands = s[,5],
                     Upgo = s[,6],
                     Armcurl = s[,7],
                     StepsNorm = s[,11],
                     StandNorm = s[,9],
                     UpgoNorm = s[,12],
                     ArmcurlNorm = s[,10])
  norm  
}
