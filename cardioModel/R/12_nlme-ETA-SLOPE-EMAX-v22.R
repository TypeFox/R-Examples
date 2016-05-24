###################################################################################
###################ETAs IN BASELINE & SLOPE/EMAX w/o BOV for PERD##################
###################################################################################

myModeling12 <- function(myData, AIC.k, delay.confidence){

  #library(nlme)

  i.int <- mean(myData$RESPONSE, na.rm = TRUE)
  i.amp <- sd(myData$RESPONSE, na.rm = TRUE)
  i.tpeak1 <- 0
  i.tpeak2 <- 0

  mask.na.for.fNTAFD <-
    !(myData[, "EXPOSURE"] %in% NA) &
    !(myData[, "RESPONSE"] %in% NA) &  !(myData[, "fNTAFD"] %in% NA)
  i.fNTAFD <- rep(i.int,length(unique(myData$fNTAFD[mask.na.for.fNTAFD])))
  i.slo <- 0
  i.emax <- i.amp
  i.ec50 <- log(mean(myData$EXPOSURE, na.rm = TRUE))
  i.hill <- log(1)

  dfEmpty <- data.frame(MODEL=NA, AIC=NA, SLOPE=NA, VAR.SLOPE=NA, EMAX=NA, VAR.EMAX=NA,
                        EC50=NA, VAR.EC50=NA, COV.EMAX.EC50=NA, HILL=NA, VAR.HILL=NA,
                        COV.EMAX.HILL=NA, COV.EC50.HILL=NA, DRUG.EFFECT.DELAY=NA,
                        ERROR.MESSAGE=NA)

  try(stop(""),TRUE)
  results.nlme.ETA.BL.SLO.EMAX.noBOV <- data.frame()
  saved.models = list()
  
  #################Estimated E0############################################

  ##Mod02-E0+Emax
  tmp <- dfEmpty
  tmp$MODEL <- "31, model = e0 + emax, bsv = e0 + emax"
  suppressWarnings(try({
    my.nlme <- nlme(RESPONSE~.my.intercept.fun(Int)+.my.Emax.fun(Emax, EC50, Conc=EXPOSURE),
                    data=myData,
                    fixed=Int+Emax+EC50~1,
                    random=list(ID=Int+Emax~1),
                    start=c(Int=i.int, Emax=i.emax, EC50=i.ec50),
                    na.action=na.omit)
    extract <- extract.param.cov(my.nlme, c("Emax", "EC50"))
    extract.num <- as.numeric(extract)
    names(extract.num) <- names(extract)
    tmp$AIC <- AIC(my.nlme, k = AIC.k)
    tmp$EMAX <- extract.num["value.Emax"]
    tmp$EC50 <- extract.num["value.EC50"]
    tmp$VAR.EMAX <- extract.num["cov.Emax.Emax"]
    tmp$VAR.EC50 <- extract.num["cov.EC50.EC50"]
    tmp$COV.EMAX.EC50 <- extract.num["cov.Emax.EC50"]
    tmp$DRUG.EFFECT.DELAY <- cder.regression(myData, my.nlme, delay.confidence)
    saved.models[[tmp$MODEL]] = my.nlme
  }, silent = T))
  tmp$ERROR.MESSAGE <- geterrmessage()
  try(stop(""),TRUE)
  results.nlme.ETA.BL.SLO.EMAX.noBOV <- rbind(results.nlme.ETA.BL.SLO.EMAX.noBOV, tmp)

  ##Mod03-E0+SigEma
  tmp <- dfEmpty
  tmp$MODEL <- "32, model = e0 + sigmoidal emax, bsv = e0 + emax"
  suppressWarnings(try({
    my.nlme <- nlme(RESPONSE~.my.intercept.fun(Int)+.my.sigEmax.fun(Emax, EC50, Hill, Conc=EXPOSURE),
                    data=myData,
                    fixed=Int+Emax+EC50+Hill~1,
                    random=list(ID=Int+Emax~1),
                    start=c(Int=i.int, Emax=i.emax, EC50=i.ec50, Hill=i.hill),
                    na.action=na.omit)
    extract <- extract.param.cov(my.nlme, c("Emax", "EC50", "Hill"))
    extract.num <- as.numeric(extract)
    names(extract.num) <- names(extract)
    tmp$AIC <- AIC(my.nlme, k = AIC.k)
    tmp$EMAX <- extract.num["value.Emax"]
    tmp$EC50 <- extract.num["value.EC50"]
    tmp$HILL <- extract.num["value.Hill"]
    tmp$VAR.EMAX <- extract.num["cov.Emax.Emax"]
    tmp$VAR.EC50 <- extract.num["cov.EC50.EC50"]
    tmp$VAR.HILL <- extract.num["cov.Hill.Hill"]
    tmp$COV.EMAX.EC50 <- extract.num["cov.Emax.EC50"]
    tmp$COV.EMAX.HILL <- extract.num["cov.Emax.Hill"]
    tmp$COV.EC50.HILL <- extract.num["cov.EC50.Hill"]
    tmp$DRUG.EFFECT.DELAY <- cder.regression(myData, my.nlme, delay.confidence)
    saved.models[[tmp$MODEL]] = my.nlme
  }, silent = T))
  tmp$ERROR.MESSAGE <- geterrmessage()
  try(stop(""),TRUE)
  results.nlme.ETA.BL.SLO.EMAX.noBOV <- rbind(results.nlme.ETA.BL.SLO.EMAX.noBOV, tmp)

  #################Estimated E0+fNTAFD##########################################

  ##Mod05-E0/fNTAFD+Emax
  tmp <- dfEmpty
  tmp$MODEL <- "33, model = e0~fNTAFD + emax, bsv = e0 + emax"
  suppressWarnings(try({
    my.nlme <- nlme(RESPONSE~.my.intercept.fun(Int)+.my.Emax.fun(Emax, EC50, Conc=EXPOSURE),
                    data=myData,
                    fixed=list(Int~fNTAFD-1, Emax+EC50~1),
                    random=list(ID=Int+Emax~1),
                    start=c(i.fNTAFD, Emax=i.emax, EC50=i.ec50),
                    na.action=na.omit)
    extract <- extract.param.cov(my.nlme, c("Emax", "EC50"))
    extract.num <- as.numeric(extract)
    names(extract.num) <- names(extract)
    tmp$AIC <- AIC(my.nlme, k = AIC.k)
    tmp$EMAX <- extract.num["value.Emax"]
    tmp$EC50 <- extract.num["value.EC50"]
    tmp$VAR.EMAX <- extract.num["cov.Emax.Emax"]
    tmp$VAR.EC50 <- extract.num["cov.EC50.EC50"]
    tmp$COV.EMAX.EC50 <- extract.num["cov.Emax.EC50"]
    tmp$DRUG.EFFECT.DELAY <- cder.regression(myData, my.nlme, delay.confidence)
    saved.models[[tmp$MODEL]] = my.nlme
  }, silent = T))
  tmp$ERROR.MESSAGE <- geterrmessage()
  try(stop(""),TRUE)
  results.nlme.ETA.BL.SLO.EMAX.noBOV <- rbind(results.nlme.ETA.BL.SLO.EMAX.noBOV, tmp)

  ##Mod06-E0/fNTAFD+SigEmax
  tmp <- dfEmpty
  tmp$MODEL <- "34, model = e0~fNTAFD + sigmoidal emax, bsv = e0 + emax"
  suppressWarnings(try({
    my.nlme <- nlme(RESPONSE~.my.intercept.fun(Int)+.my.sigEmax.fun(Emax, EC50, Hill, Conc=EXPOSURE),
                    data=myData,
                    fixed=list(Int~fNTAFD-1, Emax+EC50+Hill~1),
                    random=list(ID=Int+Emax~1),
                    start=c(i.fNTAFD, Emax=i.emax, EC50=i.ec50, Hill=i.hill),
                    na.action=na.omit)
    extract <- extract.param.cov(my.nlme, c("Emax", "EC50", "Hill"))
    extract.num <- as.numeric(extract)
    names(extract.num) <- names(extract)
    tmp$AIC <- AIC(my.nlme, k = AIC.k)
    tmp$EMAX <- extract.num["value.Emax"]
    tmp$EC50 <- extract.num["value.EC50"]
    tmp$HILL <- extract.num["value.Hill"]
    tmp$VAR.EMAX <- extract.num["cov.Emax.Emax"]
    tmp$VAR.EC50 <- extract.num["cov.EC50.EC50"]
    tmp$VAR.HILL <- extract.num["cov.Hill.Hill"]
    tmp$COV.EMAX.EC50 <- extract.num["cov.Emax.EC50"]
    tmp$COV.EMAX.HILL <- extract.num["cov.Emax.Hill"]
    tmp$COV.EC50.HILL <- extract.num["cov.EC50.Hill"]
    tmp$DRUG.EFFECT.DELAY <- cder.regression(myData, my.nlme, delay.confidence)
    saved.models[[tmp$MODEL]] = my.nlme
  }, silent = T))
  tmp$ERROR.MESSAGE <- geterrmessage()
  try(stop(""),TRUE)
  results.nlme.ETA.BL.SLO.EMAX.noBOV <- rbind(results.nlme.ETA.BL.SLO.EMAX.noBOV, tmp)

  #################Cosine Period 24h############################################

  ##Mod08-Cos24+Emax
  tmp <- dfEmpty
  tmp$MODEL <- "35, model = cosine 24 h period + emax, bsv = mean + emax"
  suppressWarnings(try({
    my.nlme <- nlme(RESPONSE~.my.cosine24.fun(Mean, Amp1, Tpeak1, TOD=TOD)+
                      .my.Emax.fun(Emax, EC50, Conc=EXPOSURE),
                    data=myData,
                    fixed=Mean+Amp1+Tpeak1+Emax+EC50~1,
                    random=list(ID=Mean+Emax~1),
                    start=c(Mean=i.int, Amp1=i.amp, Tpeak1=i.tpeak1, Emax=i.emax,
                            EC50=i.ec50),
                    na.action=na.omit)
    extract <- extract.param.cov(my.nlme, c("Emax", "EC50"))
    extract.num <- as.numeric(extract)
    names(extract.num) <- names(extract)
    tmp$AIC <- AIC(my.nlme, k = AIC.k)
    tmp$EMAX <- extract.num["value.Emax"]
    tmp$EC50 <- extract.num["value.EC50"]
    tmp$VAR.EMAX <- extract.num["cov.Emax.Emax"]
    tmp$VAR.EC50 <- extract.num["cov.EC50.EC50"]
    tmp$COV.EMAX.EC50 <- extract.num["cov.Emax.EC50"]
    tmp$DRUG.EFFECT.DELAY <- cder.regression(myData, my.nlme, delay.confidence)
    saved.models[[tmp$MODEL]] = my.nlme
  }, silent = T))
  tmp$ERROR.MESSAGE <- geterrmessage()
  try(stop(""),TRUE)
  results.nlme.ETA.BL.SLO.EMAX.noBOV <- rbind(results.nlme.ETA.BL.SLO.EMAX.noBOV, tmp)

  ##Mod09-Cos24+SigEmax
  tmp <- dfEmpty
  tmp$MODEL <- "36, model = cosine 24 h period + sigmoidal emax, bsv = mean + emax"
  suppressWarnings(try({
    my.nlme <- nlme(RESPONSE~.my.cosine24.fun(Mean, Amp1, Tpeak1, TOD=TOD)+
                      .my.sigEmax.fun(Emax, EC50, Hill, Conc=EXPOSURE),
                    data=myData,
                    fixed=Mean+Amp1+Tpeak1+Emax+EC50+Hill~1,
                    random=list(ID=Mean+Emax~1),
                    start=c(Mean=i.int, Amp1=i.amp, Tpeak1=i.tpeak1, Emax=i.emax,
                            EC50=i.ec50, Hill=i.hill),
                    na.action=na.omit)
    extract <- extract.param.cov(my.nlme, c("Emax", "EC50", "Hill"))
    extract.num <- as.numeric(extract)
    names(extract.num) <- names(extract)
    tmp$AIC <- AIC(my.nlme, k = AIC.k)
    tmp$EMAX <- extract.num["value.Emax"]
    tmp$EC50 <- extract.num["value.EC50"]
    tmp$HILL <- extract.num["value.Hill"]
    tmp$VAR.EMAX <- extract.num["cov.Emax.Emax"]
    tmp$VAR.EC50 <- extract.num["cov.EC50.EC50"]
    tmp$VAR.HILL <- extract.num["cov.Hill.Hill"]
    tmp$COV.EMAX.EC50 <- extract.num["cov.Emax.EC50"]
    tmp$COV.EMAX.HILL <- extract.num["cov.Emax.Hill"]
    tmp$COV.EC50.HILL <- extract.num["cov.EC50.Hill"]
    tmp$DRUG.EFFECT.DELAY <- cder.regression(myData, my.nlme, delay.confidence)
    saved.models[[tmp$MODEL]] = my.nlme
  }, silent = T))
  tmp$ERROR.MESSAGE <- geterrmessage()
  try(stop(""),TRUE)
  results.nlme.ETA.BL.SLO.EMAX.noBOV <- rbind(results.nlme.ETA.BL.SLO.EMAX.noBOV, tmp)

  #################Cosine Period 12h############################################

  ##Mod11-Cos12+Emax
  tmp <- dfEmpty
  tmp$MODEL <- "37, model = cosine 12 h period + emax, bsv = mean + emax"
  suppressWarnings(try({
    my.nlme <- nlme(RESPONSE~.my.cosine12.fun(Mean, Amp2, Tpeak2, TOD=TOD)+
                      .my.Emax.fun(Emax, EC50, Conc=EXPOSURE),
                    data=myData,
                    fixed=Mean+Amp2+Tpeak2+Emax+EC50~1,
                    random=list(ID=Mean+Emax~1),
                    start=c(Mean=i.int, Amp2=i.amp, Tpeak2=i.tpeak2, Emax=i.emax, EC50=i.ec50),
                    na.action=na.omit)
    extract <- extract.param.cov(my.nlme, c("Emax", "EC50"))
    extract.num <- as.numeric(extract)
    names(extract.num) <- names(extract)
    tmp$AIC <- AIC(my.nlme, k = AIC.k)
    tmp$EMAX <- extract.num["value.Emax"]
    tmp$EC50 <- extract.num["value.EC50"]
    tmp$VAR.EMAX <- extract.num["cov.Emax.Emax"]
    tmp$VAR.EC50 <- extract.num["cov.EC50.EC50"]
    tmp$COV.EMAX.EC50 <- extract.num["cov.Emax.EC50"]
    tmp$DRUG.EFFECT.DELAY <- cder.regression(myData, my.nlme, delay.confidence)
    saved.models[[tmp$MODEL]] = my.nlme
  }, silent = T))
  tmp$ERROR.MESSAGE <- geterrmessage()
  try(stop(""),TRUE)
  results.nlme.ETA.BL.SLO.EMAX.noBOV <- rbind(results.nlme.ETA.BL.SLO.EMAX.noBOV, tmp)

  ##Mod12-Cos12+SigEmax
  tmp <- dfEmpty
  tmp$MODEL <- "38, model = cosine 12 h period + sigmoidal emax, bsv = mean + emax"
  suppressWarnings(try({
    my.nlme <- nlme(RESPONSE~.my.cosine12.fun(Mean, Amp2, Tpeak2, TOD=TOD)+
                      .my.sigEmax.fun(Emax, EC50, Hill, Conc=EXPOSURE),
                    data=myData,
                    fixed=Mean+Amp2+Tpeak2+Emax+EC50+Hill~1,
                    random=list(ID=Mean+Emax~1),
                    start=c(Mean=i.int, Amp2=i.amp, Tpeak2=i.tpeak2, Emax=i.emax,
                            EC50=i.ec50, Hill=i.hill),
                    na.action=na.omit)
    extract <- extract.param.cov(my.nlme, c("Emax", "EC50", "Hill"))
    extract.num <- as.numeric(extract)
    names(extract.num) <- names(extract)
    tmp$AIC <- AIC(my.nlme, k = AIC.k)
    tmp$EMAX <- extract.num["value.Emax"]
    tmp$EC50 <- extract.num["value.EC50"]
    tmp$HILL <- extract.num["value.Hill"]
    tmp$VAR.EMAX <- extract.num["cov.Emax.Emax"]
    tmp$VAR.EC50 <- extract.num["cov.EC50.EC50"]
    tmp$VAR.HILL <- extract.num["cov.Hill.Hill"]
    tmp$COV.EMAX.EC50 <- extract.num["cov.Emax.EC50"]
    tmp$COV.EMAX.HILL <- extract.num["cov.Emax.Hill"]
    tmp$COV.EC50.HILL <- extract.num["cov.EC50.Hill"]
    tmp$DRUG.EFFECT.DELAY <- cder.regression(myData, my.nlme, delay.confidence)
    saved.models[[tmp$MODEL]] = my.nlme
  }, silent = T))
  tmp$ERROR.MESSAGE <- geterrmessage()
  try(stop(""),TRUE)
  results.nlme.ETA.BL.SLO.EMAX.noBOV <- rbind(results.nlme.ETA.BL.SLO.EMAX.noBOV, tmp)

  #################Double Cosine############################################

  ##Mod14-Dcos+Emax
  tmp <- dfEmpty
  tmp$MODEL <- "39, model = double cosine + emax, bsv = mean + emax"
  suppressWarnings(try({
    my.nlme <- nlme(RESPONSE~.my.dcosine.fun(Mean, Amp1, Tpeak1, Amp2, Tpeak2, TOD=TOD)+
                      .my.Emax.fun(Emax, EC50, Conc=EXPOSURE),
                    data=myData,
                    fixed=Mean+Amp1+Tpeak1+Amp2+Tpeak2+Emax+EC50~1,
                    random=list(ID=Mean+Emax~1),
                    start=c(Mean=i.int, Amp1=i.amp, Tpeak1=i.tpeak1, Amp2=i.amp,
                            Tpeak2=i.tpeak2, Emax=i.emax, EC50=i.ec50),
                    na.action=na.omit)
    extract <- extract.param.cov(my.nlme, c("Emax", "EC50"))
    extract.num <- as.numeric(extract)
    names(extract.num) <- names(extract)
    tmp$AIC <- AIC(my.nlme, k = AIC.k)
    tmp$EMAX <- extract.num["value.Emax"]
    tmp$EC50 <- extract.num["value.EC50"]
    tmp$VAR.EMAX <- extract.num["cov.Emax.Emax"]
    tmp$VAR.EC50 <- extract.num["cov.EC50.EC50"]
    tmp$COV.EMAX.EC50 <- extract.num["cov.Emax.EC50"]
    tmp$DRUG.EFFECT.DELAY <- cder.regression(myData, my.nlme, delay.confidence)
    saved.models[[tmp$MODEL]] = my.nlme
  }, silent = T))
  tmp$ERROR.MESSAGE <- geterrmessage()
  try(stop(""),TRUE)
  results.nlme.ETA.BL.SLO.EMAX.noBOV <- rbind(results.nlme.ETA.BL.SLO.EMAX.noBOV, tmp)

  ##Mod15-Dcos+SigEmax
  tmp <- dfEmpty
  tmp$MODEL <- "40, model = double cosine + sigmoidal emax, bsv = mean + emax"
  suppressWarnings(try({
    my.nlme <- nlme(RESPONSE~.my.dcosine.fun(Mean, Amp1, Tpeak1, Amp2, Tpeak2, TOD=TOD)+
                      .my.sigEmax.fun(Emax, EC50, Hill, Conc=EXPOSURE),
                    data=myData,
                    fixed=Mean+Amp1+Tpeak1+Amp2+Tpeak2+Emax+EC50+Hill~1,
                    random=list(ID=Mean+Emax~1),
                    start=c(Mean=i.int, Amp1=i.amp, Tpeak1=i.tpeak1, Amp2=i.amp,
                            Tpeak2=i.tpeak2, Emax=i.emax, EC50=i.ec50, Hill=i.hill),
                    na.action=na.omit)
    extract <- extract.param.cov(my.nlme, c("Emax", "EC50", "Hill"))
    extract.num <- as.numeric(extract)
    names(extract.num) <- names(extract)
    tmp$AIC <- AIC(my.nlme, k = AIC.k)
    tmp$EMAX <- extract.num["value.Emax"]
    tmp$EC50 <- extract.num["value.EC50"]
    tmp$HILL <- extract.num["value.Hill"]
    tmp$VAR.EMAX <- extract.num["cov.Emax.Emax"]
    tmp$VAR.EC50 <- extract.num["cov.EC50.EC50"]
    tmp$VAR.HILL <- extract.num["cov.Hill.Hill"]
    tmp$COV.EMAX.EC50 <- extract.num["cov.Emax.EC50"]
    tmp$COV.EMAX.HILL <- extract.num["cov.Emax.Hill"]
    tmp$COV.EC50.HILL <- extract.num["cov.EC50.Hill"]
    tmp$DRUG.EFFECT.DELAY <- cder.regression(myData, my.nlme, delay.confidence)
    saved.models[[tmp$MODEL]] = my.nlme
  }, silent = T))
  tmp$ERROR.MESSAGE <- geterrmessage()
  try(stop(""),TRUE)
  results.nlme.ETA.BL.SLO.EMAX.noBOV <- rbind(results.nlme.ETA.BL.SLO.EMAX.noBOV, tmp)

  #Replace fake error by none
  mask.fake.error <- grep("try", results.nlme.ETA.BL.SLO.EMAX.noBOV$ERROR.MESSAGE)
  results.nlme.ETA.BL.SLO.EMAX.noBOV$ERROR.MESSAGE[mask.fake.error] <-
    "none"

  #Return results
  return(list(results.nlme.ETA.BL.SLO.EMAX.noBOV, saved.models))
}
