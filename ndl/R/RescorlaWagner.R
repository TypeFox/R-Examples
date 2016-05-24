RescorlaWagner = function(
  cuesOutcomes, 
  traceCue = "h",
  traceOutcome = "hand",
  nruns=1, 
  random = TRUE,
  randomOrder = NA,
  alpha=0.1, lambda=1, beta1=0.1,beta2=0.1) {

  cues =  unique(unlist(strsplit(as.character(cuesOutcomes$Cues), "_")))
  outcomes = unique(unlist(strsplit(as.character(cuesOutcomes$Outcomes), "_")))

  weightvec = rep(0, length(cues))
  names(weightvec) = cues

  res = vector(mode="numeric")
  # ultimately, res will have nruns * sum(cuesOutcomes$Frequency) elements

  dfr = data.frame(
    cuevector = rep(cuesOutcomes$Cues, cuesOutcomes$Frequency),
    outcomevector = rep(cuesOutcomes$Outcomes, cuesOutcomes$Frequency),
    stringsAsFactors=FALSE
  )

  theOrder = NA
  cnt = 0
  for (run in 1:nruns) {
    if (random) {
      if (is.na(randomOrder[1])) {
        theOrder = sample(1:nrow(dfr))
      } else {
        theOrder = randomOrder
      }
      dfr = dfr[theOrder,]
    }
    for (i in 1:nrow(dfr)) {
      currentCues = unlist(strsplit(dfr$cuevector[i],"_"))
      currentOutcomes = unlist(strsplit(dfr$outcomevector[i],"_"))
      Vtotal = sum(weightvec[currentCues])
      if (traceOutcome %in% currentOutcomes) {
        Lambda = lambda
      } else {
        Lambda = 0
      }
      #for (j in 1:length(weightvec)) {  
      #  if (names(weightvec)[j] %in% currentCues) {
      #    weightvec[j] = weightvec[j] + alpha*beta1*(Lambda-Vtotal)
      #  }
      #}
      weightvec[currentCues] = weightvec[currentCues]+alpha*beta1*(Lambda-Vtotal)
      cnt = cnt + 1 
      res[cnt] = weightvec[traceCue]
    }
  }

  equilibriumWeights = estimateWeights(cuesOutcomes)
  eqw = equilibriumWeights[traceCue, traceOutcome]

  result <- (list(weightvector=res, equilibriumWeight=eqw,
          traceCue=traceCue, traceOutcome=traceOutcome, randomOrder=theOrder))
  class(result) <- "RescorlaWagner"

  return(result)
}
