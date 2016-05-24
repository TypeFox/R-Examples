###################################################################################
###################ETAs IN BASELINE w/o BOV for PERD###############################
###################################################################################

myModeling18 <- function(myData, AIC.k, delay.confidence){

  i.int <- mean(myData$RESPONSE, na.rm = TRUE)
  i.amp <- sd(myData$RESPONSE, na.rm = TRUE)
  i.tpeak1 <- 0
  i.tpeak2 <- 0

  mask.na.for.fNTAFD <-
    !(myData[, "EXPOSURE"] %in% NA) &
    !(myData[, "RESPONSE"] %in% NA) &  !(myData[, "fNTAFD"] %in% NA)
  i.fNTAFD <- rep(i.int,length(unique(myData$fNTAFD[mask.na.for.fNTAFD])))

  dfEmpty <- data.frame(MODEL=NA, AIC=NA, SLOPE=NA, VAR.SLOPE=NA, EMAX=NA, VAR.EMAX=NA,
                        EC50=NA, VAR.EC50=NA, COV.EMAX.EC50=NA, HILL=NA, VAR.HILL=NA,
                        COV.EMAX.HILL=NA, COV.EC50.HILL=NA, DRUG.EFFECT.DELAY=NA,
                        ERROR.MESSAGE=NA)

  try(stop(""),TRUE)
  results.nlme.ETA.BL.noBOV <- data.frame()
  saved.models = list()
  
  #################Estimated E0############################################

  # Mod01-E0+Slope
  tmp <- dfEmpty
  tmp$MODEL <- "101, model = e0, bsv = e0"
  suppressWarnings(try({
    my.nlme <- nlme(RESPONSE~.my.intercept.fun(Int),
                    data=myData,
                    fixed=Int~1,
                    random=list(ID=Int~1),
                    start=c(Int=i.int),
                    na.action=na.omit)
    tmp$AIC <- AIC(my.nlme, k = AIC.k)
    saved.models[[tmp$MODEL]] = my.nlme
  }, silent = T))
  tmp$ERROR.MESSAGE <- geterrmessage()
  try(stop(""),TRUE)
  results.nlme.ETA.BL.noBOV <- rbind(results.nlme.ETA.BL.noBOV, tmp)

  #################Estimated E0+fNTAFD##########################################

  #Mod04-E0/fNTAFD+Slope
  tmp <- dfEmpty
  tmp$MODEL <- "102, model = e0~fNTAFD, bsv = e0"
  suppressWarnings(try({
    my.nlme <- nlme(RESPONSE~.my.intercept.fun(Int),
                    data=myData,
                    fixed=list(Int~fNTAFD-1),
                    random=list(ID=Int~1),
                    start=c(i.fNTAFD),
                    na.action=na.omit)
    tmp$AIC <- AIC(my.nlme, k = AIC.k)
    saved.models[[tmp$MODEL]] = my.nlme
  }, silent = T))
  tmp$ERROR.MESSAGE <- geterrmessage()
  try(stop(""),TRUE)
  results.nlme.ETA.BL.noBOV <- rbind(results.nlme.ETA.BL.noBOV, tmp)

  #################Cosine Period 24h############################################

  ##Mod07-Cos24+Slope
  tmp <- dfEmpty
  tmp$MODEL <- "103, model = cosine 24 h period, bsv = mean"
  suppressWarnings(try({
    my.nlme <- nlme(RESPONSE~.my.cosine24.fun(Mean, Amp1, Tpeak1, TOD=TOD),
                    data=myData,
                    fixed=Mean+Amp1+Tpeak1~1,
                    random=list(ID=Mean~1),
                    start=c(Mean=i.int, Amp1=i.amp, Tpeak1=i.tpeak1),
                    na.action=na.omit)
    tmp$AIC <- AIC(my.nlme, k = AIC.k)
    saved.models[[tmp$MODEL]] = my.nlme
  }, silent = T))
  tmp$ERROR.MESSAGE <- geterrmessage()
  try(stop(""),TRUE)
  results.nlme.ETA.BL.noBOV <- rbind(results.nlme.ETA.BL.noBOV, tmp)


  #################Cosine Period 12h############################################

  ##Mod10-Cos12+Slope
  tmp <- dfEmpty
  tmp$MODEL <- "104, model = cosine 12 h period, bsv = mean"
  suppressWarnings(try({
    my.nlme <- nlme(RESPONSE~.my.cosine12.fun(Mean, Amp2, Tpeak2, TOD=TOD),
                    data=myData,
                    fixed=Mean+Amp2+Tpeak2~1,
                    random=list(ID=Mean~1),
                    start=c(Mean=i.int, Amp2=i.amp, Tpeak2=i.tpeak2),
                    na.action=na.omit)
    tmp$AIC <- AIC(my.nlme, k = AIC.k)
    saved.models[[tmp$MODEL]] = my.nlme
  }, silent = T))
  tmp$ERROR.MESSAGE <- geterrmessage()
  try(stop(""),TRUE)
  results.nlme.ETA.BL.noBOV <- rbind(results.nlme.ETA.BL.noBOV, tmp)


  #################Double Cosine############################################

  ##Mod13-Dcos+Slope
  tmp <- dfEmpty
  tmp$MODEL <- "105, model = double cosine, bsv = mean"
  suppressWarnings(try({
    my.nlme <- nlme(RESPONSE~.my.dcosine.fun(Mean, Amp1, Tpeak1, Amp2, Tpeak2, TOD=TOD),
                    data=myData,
                    fixed=Mean+Amp1+Tpeak1+Amp2+Tpeak2~1,
                    random=list(ID=Mean~1),
                    start=c(Mean=i.int, Amp1=i.amp, Tpeak1=i.tpeak1, Amp2=i.amp,
                            Tpeak2=i.tpeak2),
                    na.action=na.omit)
    tmp$AIC <- AIC(my.nlme, k = AIC.k)
    saved.models[[tmp$MODEL]] = my.nlme
  }, silent = T))
  tmp$ERROR.MESSAGE <- geterrmessage()
  try(stop(""),TRUE)
  results.nlme.ETA.BL.noBOV <- rbind(results.nlme.ETA.BL.noBOV, tmp)

  #Replace fake error by none
  mask.fake.error <- grep("try", results.nlme.ETA.BL.noBOV$ERROR.MESSAGE)
  results.nlme.ETA.BL.noBOV$ERROR.MESSAGE[mask.fake.error] <-
    "none"

  #Return results
  return(list(results.nlme.ETA.BL.noBOV, saved.models))
}
