##' testthat test to check multiple varying parameters
##' on all of the functions of wrassp somehow breaks the
##' functions 
##'
##' @author Raphael Winkelmann
context("test signal processing functions")

###################################
# acfana
test_that("acfana doesn't break due to varying parameters", {

  nrOfRandomCalls = 10

  wavFiles <- list.files(system.file("extdata", package = "wrassp"), pattern = glob2rx("*.wav"), full.names = TRUE)

  posValsBeginTime = list(0, 0.0001, 0.1, 0.5)
  posValsCenterTime = list(TRUE, FALSE)
  posValsEndTime = list(0, 0.7, 0.700001, 1, 1.2)
  posValsWindowShift = list(5) #list(1, 5, 10, 7)
  posValsWindowSize = list(20) #???
  posValsEffectiveLength = list(TRUE) #list(TRUE, FALSE)
  posValsWindow = list("BLACKMAN")#as.list(AsspWindowTypes())
  posValsAnalysisOrder = list(0) #???
  posValsEnergyNormalization = list(FALSE) #list(TRUE, FALSE)
  posValsLengthNormalization = list(FALSE) #list(TRUE, FALSE)

  for(i in 1:nrOfRandomCalls){
    params = list(listOfFiles=sample(wavFiles, 1)[[1]], optLogFilePath=NULL, 
                  beginTime=sample(posValsBeginTime, 1)[[1]], centerTime=sample(posValsCenterTime, 1)[[1]],
                  endTime=sample(posValsEndTime, 1)[[1]], windowShift=sample(posValsWindowShift, 1)[[1]], 
                  windowSize=sample(posValsWindowSize, 1)[[1]], effectiveLength=sample(posValsEffectiveLength, 1)[[1]], 
                  window=sample(posValsWindow, 1)[[1]], analysisOrder=sample(posValsAnalysisOrder, 1)[[1]], 
                  energyNormalization=sample(posValsEnergyNormalization, 1)[[1]], lengthNormalization=sample(posValsLengthNormalization, 1)[[1]], 
                  explicitExt=NULL, outputDirectory=NULL,
                  toFile=FALSE, forceToLog=useWrasspLogger)
    # print(params)
    res = do.call(acfana,as.list(params))

    expect_that(class(res), equals("AsspDataObj"))
  }

})

##################################
# afdiff
test_that("afdiff doesn't break due to varying parameters", {

  nrOfRandomCalls = 10

  wavFiles <- list.files(system.file("extdata", package = "wrassp"), pattern = glob2rx("*.wav"), full.names = TRUE)

  posValsComputeBackwardDifference = list(FALSE) #list(TRUE, FALSE)
  posValsComputeCentralDifference = list(FALSE) #list(TRUE, FALSE)
  posValsChannel = list(1)

  for(i in 1:nrOfRandomCalls){
    params = list(listOfFiles=sample(wavFiles, 1)[[1]], optLogFilePath=NULL, 
                  computeBackwardDifferenc=sample(posValsComputeBackwardDifference,1)[[1]], computeCentralDifference=sample(posValsComputeCentralDifference,1)[[1]],
                  channel=sample(posValsChannel,1)[[1]], toFile=FALSE,
                  explicitExt=NULL, outputDirectory=NULL,
                  forceToLog=useWrasspLogger)
    # print(params)
    res = do.call(afdiff,as.list(params))
    expect_that(class(res), equals("AsspDataObj"))
  }

})

##################################
# affilter
test_that("affilter doesn't break due to varying parameters", {

  nrOfRandomCalls = 10

  wavFiles <- list.files(system.file("extdata", package = "wrassp"), pattern = glob2rx("*.wav"), full.names = TRUE)

  posValsHighPass=list(4000)
  posValsLowPass=list(0)
  posValsStopBand=list(96)
  posValsTransition=list(250)
  posValsUseIIR=list(FALSE)
  posValsNumIIRsections=list(4)

  for(i in 1:nrOfRandomCalls){
    params = list(listOfFiles=sample(wavFiles, 1)[[1]], optLogFilePath=NULL,
                  highPass=sample(posValsHighPass,1)[[1]], lowPass=sample(posValsLowPass,1)[[1]], 
                  stopBand=sample(posValsStopBand,1)[[1]], transition=sample(posValsTransition,1)[[1]], 
                  useIIR=sample(posValsUseIIR,1)[[1]], numIIRsections=sample(posValsNumIIRsections,1)[[1]],
                  toFile=FALSE, explicitExt=NULL,
                  outputDirectory=NULL, forceToLog=useWrasspLogger)
  
    # print(params)
    res = do.call(affilter, as.list(params))
    expect_that(class(res), equals("AsspDataObj"))
  }
})

##################################
# cepstrum
test_that("cepstrum doesn't break due to varying parameters", {

  nrOfRandomCalls = 10

  wavFiles <- list.files(system.file("extdata", package = "wrassp"), pattern = glob2rx("*.wav"), full.names = TRUE)

  posValsBeginTime=list(0, 0.0001, 0.1, 0.5)
  posValsCenterTime=list(TRUE, FALSE)
  posValsEndTime=list(0, 0.7, 0.700001, 1, 1.2)
  posValsResolution=list(40)
  posValsFftLength=list(0)
  posValsWindowShift=list(5)
  posValsWindow=list("BLACKMAN")#as.list(AsspWindowTypes())

  for(i in 1:nrOfRandomCalls){
    params = list(listOfFiles=sample(wavFiles, 1)[[1]], optLogFilePath=NULL, 
                  beginTime=sample(posValsBeginTime,1)[[1]], centerTime=sample(posValsCenterTime,1)[[1]],
                  endTime=sample(posValsEndTime,1)[[1]], resolution=sample(posValsResolution,1)[[1]],
                  fftLength=sample(posValsFftLength,1)[[1]], windowShift=sample(posValsWindowShift,1)[[1]],
                  window=sample(posValsWindow,1)[[1]], toFile=FALSE,
                  explicitExt=NULL, outputDirectory=NULL,
                  forceToLog=useWrasspLogger)
    
    # print(params)
    res = do.call(cepstrum, as.list(params))
    expect_that(class(res), equals("AsspDataObj"))
  }


})

##################################
# cssSpectrum
test_that("cssSpectrum doesn't break due to varying parameters", {

  nrOfRandomCalls = 10

  wavFiles <- list.files(system.file("extdata", package = "wrassp"), pattern = glob2rx("*.wav"), full.names = TRUE)

  posValsBeginTime=list(0, 0.0001, 0.1, 0.5)
  posValsCenterTime=list(TRUE, FALSE)
  posValsEndTime=list(0, 0.7, 0.700001, 1, 1.2)
  posValsResolution=list(40)
  posValsFftLength=list(0)
  posValsWindowShift=list(5)
  posValsWindow=list("BLACKMAN") #as.list(AsspWindowTypes())
  posValsNumCeps=list(0)

  for (i in 1:nrOfRandomCalls) {
    params = list(listOfFiles=sample(wavFiles, 1)[[1]], optLogFilePath=NULL, 
                  beginTime=sample(posValsBeginTime,1)[[1]], centerTime=sample(posValsCenterTime,1)[[1]], 
                  endTime=sample(posValsEndTime,1)[[1]], resolution=sample(posValsResolution,1)[[1]], 
                  fftLength=sample(posValsFftLength,1)[[1]], windowShift=sample(posValsWindowShift,1)[[1]], 
                  window=sample(posValsWindow,1)[[1]], numCeps=sample(posValsNumCeps,1)[[1]], 
                  toFile=FALSE, explicitExt=NULL, 
                  outputDirectory=NULL, forceToLog=useWrasspLogger)
        # print(params)
    res = do.call(cssSpectrum, as.list(params))
    expect_that(class(res), equals("AsspDataObj"))
  }
})

##################################
# dftSpectrum
test_that("dftSpectrum doesn't break due to varying parameters", {

  nrOfRandomCalls = 10

  wavFiles <- list.files(system.file("extdata", package = "wrassp"), pattern = glob2rx("*.wav"), full.names = TRUE)

  posValsBeginTime=list(0, 0.0001, 0.1, 0.5)
  posValsCenterTime=list(TRUE, FALSE)
  posValsEndTime=list(0, 0.7, 0.700001, 1, 1.2)
  posValsResolution=list(40)
  posValsFftLength=list(0)
  posValsWindowShift=list(5)
  posValsWindow=list("BLACKMAN") #as.list(AsspWindowTypes())
  posValsBandwidth=list(0)

  for (i in 1:nrOfRandomCalls) {
    params = list(listOfFiles=sample(wavFiles, 1)[[1]], optLogFilePath=NULL, 
                  beginTime=posValsBeginTime[[1]], centerTime=posValsCenterTime[[1]], 
                  endTime=posValsEndTime[[1]], resolution=posValsResolution[[1]], 
                  fftLength=posValsFftLength[[1]], windowShift=posValsWindowShift[[1]], 
                  window=posValsWindow[[1]], bandwidth=posValsBandwidth[[1]], 
                  toFile=FALSE, explicitExt=NULL, 
                  outputDirectory=NULL, forceToLog=useWrasspLogger)
  
    # print(params)
    res = do.call(dftSpectrum, as.list(params))
    expect_that(class(res), equals("AsspDataObj"))

  }  
})

##################################
# forest
test_that("forest doesn't break due to varying parameters", {

  nrOfRandomCalls = 10

  wavFiles <- list.files(system.file("extdata", package = "wrassp"), pattern = glob2rx("*.wav"), full.names = TRUE)

  posValsBeginTime=list(0, 0.0001, 0.1, 0.5)
  posValsEndTime=list(0, 0.7, 0.700001, 1, 1.2)
  posValsWindowShift=list(5)
  posValsWindowSize=list(20)
  posValsEffectiveLength=list(TRUE)
  posValsNominalF1=list(500)
  posValsGender=list("m")
  posValsEstimate=list(FALSE)
  posValsOrder=list(0)
  posValsIncrOrder=list(0)
  posValsNumFormants=list(4)
  posValsWindow=list("BLACKMAN")#as.list(AsspWindowTypes())
  posValsPreemphasis=list(-0.8)
  
  for (i in 1:nrOfRandomCalls) {
    params = list(listOfFiles=sample(wavFiles, 1)[[1]], optLogFilePath=NULL, 
                  beginTime=sample(posValsBeginTime,1)[[1]], endTime=sample(posValsEndTime,1)[[1]], 
                  windowShift=sample(posValsWindowShift,1)[[1]], windowSize=sample(posValsWindowSize,1)[[1]], 
                  effectiveLength=sample(posValsEffectiveLength,1)[[1]], nominalF1=sample(posValsNominalF1,1)[[1]], 
                  gender=sample(posValsGender,1)[[1]], estimate=sample(posValsEstimate,1)[[1]], 
                  order=sample(posValsOrder,1)[[1]], incrOrder=sample(posValsIncrOrder,1)[[1]], 
                  numFormants=sample(posValsNumFormants,1)[[1]], window=sample(posValsWindow,1)[[1]], 
                  preemphasis=sample(posValsPreemphasis,1)[[1]], toFile=FALSE, 
                  explicitExt=NULL, outputDirectory=NULL, 
                  forceToLog=useWrasspLogger)

    # print(params)
    res = do.call(forest, as.list(params))
    expect_that(class(res), equals("AsspDataObj"))
  }
})

##################################
# ksvF0
test_that("ksvF0 doesn't break due to varying parameters", {
  
  nrOfRandomCalls = 10

  wavFiles <- list.files(system.file("extdata", package = "wrassp"), pattern = glob2rx("*.wav"), full.names = TRUE)

  posValsBeginTime=list(0, 0.0001, 0.1, 0.5)
  posValsEndTime=list(0, 0.7, 0.700001, 1, 1.2)
  posValsWindowShift=list(5)
  posValsGender=list("u")
  posValsMaxF=list(600)
  posValsMinF=list(50)
  posValsMinAmp=list(50)
  posValsMaxZCR=list(3000)

  for (i in 1:nrOfRandomCalls) {
    params = list(listOfFiles=sample(wavFiles, 1)[[1]], optLogFilePath=NULL, 
                  beginTime=sample(posValsBeginTime,1)[[1]], endTime=sample(posValsEndTime,1)[[1]], 
                  windowShift=sample(posValsWindowShift,1)[[1]], gender=sample(posValsGender,1)[[1]], 
                  maxF=sample(posValsMaxF,1)[[1]], minF=sample(posValsMinF,1)[[1]], 
                  minAmp=sample(posValsMinAmp,1)[[1]], maxZCR=sample(posValsMaxZCR,1)[[1]], 
                  toFile=FALSE, explicitExt=NULL, 
                  outputDirectory=NULL, forceToLog=useWrasspLogger)

        # print(params)
    res = do.call(ksvF0, as.list(params))
    expect_that(class(res), equals("AsspDataObj"))
  }
})

##################################
# lpsSpectrum
test_that("lpsSpectrum doesn't break due to varying parameters", {
  
  nrOfRandomCalls = 10

  wavFiles <- list.files(system.file("extdata", package = "wrassp"), pattern = glob2rx("*.wav"), full.names = TRUE)
  
  posValsBeginTime=list(0, 0.0001, 0.1, 0.5)
  posValsCenterTime=list(TRUE, FALSE)
  posValsEndTime=list(0, 0.7, 0.700001, 1, 1.2)
  posValsResolution=list(40)
  posValsFftLength=list(0)
  posValsWindowSize=list(20)
  posValsWindowShift=list(5)
  posValsWindow=list("BLACKMAN") #as.list(AsspWindowTypes())
  posValsOrder=list(0)
  posValsPreemphasis=list(-0.95)
  posValsDeemphasize=list(TRUE)

  for (i in 1:nrOfRandomCalls) {
    params = list(listOfFiles=sample(wavFiles, 1)[[1]], optLogFilePath=NULL, 
                  beginTime=sample(posValsBeginTime,1)[[1]], centerTime=sample(posValsCenterTime,1)[[1]], 
                  endTime=sample(posValsEndTime,1)[[1]], resolution=sample(posValsResolution,1)[[1]], 
                  fftLength=sample(posValsFftLength,1)[[1]], windowSize=sample(posValsWindowSize,1)[[1]], 
                  windowShift=sample(posValsWindowShift,1)[[1]], window=sample(posValsWindow,1)[[1]], 
                  order=sample(posValsOrder,1)[[1]], preemphasis=sample(posValsPreemphasis,1)[[1]], 
                  deemphasize=sample(posValsDeemphasize,1)[[1]], toFile=FALSE, 
                  explicitExt=NULL, outputDirectory=NULL, 
                  forceToLog=useWrasspLogger)
    
    # print(params)
    res = do.call(lpsSpectrum, as.list(params))
    expect_that(class(res), equals("AsspDataObj"))
  }
})

##################################
# mhsF0
test_that("mhsF0 doesn't break due to varying parameters", {
  
  nrOfRandomCalls = 10

  wavFiles <- list.files(system.file("extdata", package = "wrassp"), pattern = glob2rx("*.wav"), full.names = TRUE)

  posValsBeginTime=list(0, 0.0001, 0.1, 0.5)
  posValsCenterTime=list(FALSE) #list(FALSE)
  posValsEndTime=list(0, 0.7, 0.700001, 1, 1.2)
  posValsWindowShift=list(5)
  posValsGender=list("f","m","u")
  posValsMaxF=list(600)
  posValsMinF=list(50)
  posValsMinAmp=list(50)
  posValsMinAC1=list(0.25)
  posValsMinRMS=list(18)
  posValsMaxZCR=list(3000)
  posValsMinProb=list(0.52)
  posValsPlainSpectrum=list(FALSE)

  for (i in 1:nrOfRandomCalls) {
    params = list(listOfFiles=sample(wavFiles, 1)[[1]], optLogFilePath=NULL, 
                  beginTime=sample(posValsBeginTime,1)[[1]], centerTime=sample(posValsCenterTime,1)[[1]], 
                  endTime=sample(posValsEndTime,1)[[1]], windowShift=sample(posValsWindowShift,1)[[1]], 
                  gender=sample(posValsGender,1)[[1]], maxF=sample(posValsMaxF,1)[[1]], 
                  minF=sample(posValsMinF,1)[[1]], minAmp=sample(posValsMinAmp,1)[[1]], 
                  minAC1=sample(posValsMinAC1,1)[[1]], minRMS=sample(posValsMinRMS,1)[[1]], 
                  maxZCR=sample(posValsMaxZCR,1)[[1]], minProb=sample(posValsMinProb,1)[[1]], 
                  plainSpectrum=sample(posValsPlainSpectrum,1)[[1]], toFile=FALSE, 
                  explicitExt=NULL, outputDirectory=NULL, 
                  forceToLog=useWrasspLogger)

    # print(params)
    res = do.call(mhsF0, as.list(params))
    expect_that(class(res), equals("AsspDataObj"))

  }
})

##################################
# rfcana
test_that("rfcana doesn't break due to varying parameters", {
  
  nrOfRandomCalls = 10

  wavFiles <- list.files(system.file("extdata", package = "wrassp"), pattern = glob2rx("*.wav"), full.names = TRUE)

  posValsBeginTime=list(0, 0.0001, 0.1, 0.5)
  posValsCenterTime=list(TRUE, FALSE)
  posValsEndTime=list(0, 0.7, 0.700001, 1, 1.2)
  posValsWindowShift=list(5)
  posValsWindowSize=list(20)
  posValsEffectiveLength=list(TRUE)
  posValsWindow=list("BLACKMAN")
  posValsOrder=list(0)
  posValsPreemphasis=list(-0.95)
  posValsLpType=list("RFC")

  for (i in 1:nrOfRandomCalls) {
    params = list(listOfFiles=sample(wavFiles, 1)[[1]], optLogFilePath=NULL, 
                  beginTime=sample(posValsBeginTime,1)[[1]], centerTime=sample(posValsCenterTime,1)[[1]], 
                  endTime=sample(posValsEndTime,1)[[1]], windowShift=sample(posValsWindowShift,1)[[1]], 
                  windowSize=sample(posValsWindowSize,1)[[1]], effectiveLength=sample(posValsEffectiveLength,1)[[1]], 
                  window=sample(posValsWindow,1)[[1]], order=sample(posValsOrder,1)[[1]], 
                  preemphasis=sample(posValsPreemphasis,1)[[1]], lpType=sample(posValsLpType,1)[[1]], 
                  toFile=FALSE, explicitExt=NULL, 
                  outputDirectory=NULL, forceToLog=useWrasspLogger)
        # print(params)
    res = do.call(rfcana, as.list(params))
    expect_that(class(res), equals("AsspDataObj"))

  }
})

#######################################
# rmsana
test_that("rmsana doesn't break due to varying parameters", {
    
  nrOfRandomCalls = 10

  wavFiles <- list.files(system.file("extdata", package = "wrassp"), pattern = glob2rx("*.wav"), full.names = TRUE)

  posValsBeginTime=list(0, 0.0001, 0.1, 0.5)
  posValsCenterTime=list(TRUE, FALSE)
  posValsEndTime=list(0, 0.7, 0.700001, 1, 1.2)
  posValsWindowShift=list(5)
  posValsWindowSize=list(20)
  posValsEffectiveLength=list(TRUE)
  posValsLinear=list(FALSE)
  posValsWindow=list("HAMMING")

  for (i in 1:nrOfRandomCalls) {
    params = list(listOfFiles=sample(wavFiles, 1)[[1]], optLogFilePath=NULL, 
                  beginTime=sample(posValsBeginTime,1)[[1]], centerTime=sample(posValsCenterTime,1)[[1]], 
                  endTime=sample(posValsEndTime,1)[[1]], windowShift=sample(posValsWindowShift,1)[[1]], 
                  windowSize=sample(posValsWindowSize,1)[[1]], effectiveLength=sample(posValsEffectiveLength,1)[[1]], 
                  linear=sample(posValsLinear,1)[[1]], window=sample(posValsWindow,1)[[1]], 
                  toFile=FALSE, explicitExt=NULL, 
                  outputDirectory=NULL, forceToLog=useWrasspLogger)
            # print(params)
    res = do.call(rmsana, as.list(params))
    expect_that(class(res), equals("AsspDataObj"))

  }
})

##################################
# zcrana
test_that("zcrana doesn't break due to varying parameters", {

  nrOfRandomCalls = 10

  wavFiles <- list.files(system.file("extdata", package = "wrassp"), pattern = glob2rx("*.wav"), full.names = TRUE)

  posValsBeginTime=list(0, 0.0001, 0.1, 0.5)
  posValsCenterTime=list(TRUE, FALSE)
  posValsEndTime=list(0, 0.7, 0.700001, 1, 1.2)
  posValsWindowShift=list(5)
  posValsWindowSize=list(25)

  for (i in 1:nrOfRandomCalls) {
    params = list(listOfFiles=sample(wavFiles, 1)[[1]], optLogFilePath=NULL, 
                  beginTime=sample(posValsBeginTime,1)[[1]], centerTime=sample(posValsCenterTime,1)[[1]], 
                  endTime=sample(posValsEndTime,1)[[1]], windowShift=sample(posValsWindowShift,1)[[1]], 
                  windowSize=sample(posValsWindowSize,1)[[1]], toFile=FALSE, 
                  explicitExt=NULL, outputDirectory=NULL, 
                  forceToLog=useWrasspLogger)
                # print(params)
    res = do.call(zcrana, as.list(params))
    expect_that(class(res), equals("AsspDataObj"))
  }
})
