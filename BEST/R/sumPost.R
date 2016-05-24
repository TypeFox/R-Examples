sumPost <-
function(paramSampleVec, credMass=0.95, compVal=NULL, ROPE=NULL) {
  # Gets summary information for a single parameter; 
  # called by summary.BEST; not exported.
  postSummary <- rep(NA, 11) 
  names(postSummary) <- c("mean","median","mode",
                         "hdiMass","hdiLow","hdiHigh",
                         "compVal","pcGTcompVal",
                         "ROPElow","ROPEhigh","pcInROPE")        
  postSummary["mean"] <- mean(paramSampleVec)
  postSummary["median"] <- median(paramSampleVec)
  mcmcDensity <- density(paramSampleVec)
  postSummary["mode"] <- mcmcDensity$x[which.max(mcmcDensity$y)]

  HDI <- hdi(paramSampleVec, credMass)
  postSummary["hdiMass"] <- credMass * 100
  postSummary["hdiLow"] <- HDI[1]
  postSummary["hdiHigh"] <- HDI[2]

  if (!is.null(compVal)) {
    postSummary["compVal"] <- compVal
    postSummary["pcGTcompVal"] <- mean(paramSampleVec > compVal) * 100
  }

  if (!is.null(ROPE)) {
    postSummary["ROPElow"] <- ROPE[1] 
    postSummary["ROPEhigh"] <- ROPE[2] 
    postSummary["pcInROPE"] <- mean(paramSampleVec > ROPE[1] & 
                                      paramSampleVec < ROPE[2]) * 100
  }
  return(postSummary)
}
