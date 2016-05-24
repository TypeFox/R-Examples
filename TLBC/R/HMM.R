## functions to train and test HMMs

trainHMM = function(labelDir, rf, names, combineStanding=FALSE) {
  
  # load ground truth labels
  labels = loadLabels(labelDir, names)
  labels = labels$behavior
  
  cat("training HMM\n")
  # get the unique states
  states = sort(unique(labels))
  states = states[!grepl("NULL", states)]
  symbols = sort(rf$classes)
  # compute the transition probabilities
  transProbs = computeTransProbs(labels)
  # compute the emission probabilities
  emissionProbs = computeEmissionProbs(rf)
  # compute the prior probabilities
  startProbs = computePriorProbs(labels)
  # make the HMM
  hmm = initHMM(states, symbols, startProbs, transProbs, emissionProbs)
  return(hmm)
}
testHMM = function(predDir, modelName, saveDir, names) {
  # load model
  hmm = loadModel(modelName, "hmm")
  if (length(names) == 0) {
    stop("no test data\n")
  }
  for (name in names) {
    # load predictions
    pred = loadPredictions(predDir, name)
  
    if (nrow(pred) == 0) {
      next
    }
    # apply HMM
    cat("applying HMM\n")
    
    post = posterior(hmm, pred$prediction)
    filtered = as.character(factor(max.col(t(post)), levels=1:length(hmm$States), 
                                   labels=hmm$States))
    filtered = viterbi(hmm, pred$prediction)
    
    # save predictions
    saveFile = file.path(saveDir, paste0(name, ".csv"))
    writePredictions(filtered, pred$timestamp, saveFile)
  }
}
computeTransProbs = function(stateSeq) {
  # get the list of unque states
  states = sort(unique(stateSeq))
  # remove NULL
  states = states[!grepl("NULL", states)]
  S = length(states)
  # set up the transition Probability matrix
  transProbs = matrix(.Machine$double.eps, nrow=S, ncol=S)
  rownames(transProbs) = states
  colnames(transProbs) = states
  # loop through state sequence and count transitions
  x = stateSeq[1]
  i = 2
  while (x == "NULL") {
    x = stateSeq[i]
    i = i + 1
  }
  nNull = 0
  for (ii in i:(length(stateSeq) - 1)) {
    if (stateSeq[ii] != "NULL") {
      if (nNull <= 4){
        transProbs[x,stateSeq[ii]] = transProbs[x, stateSeq[ii]] + 1
        x = stateSeq[ii]
        nNull = 0
      } else {
        nNull = 0
      }
    } else {
      nNull = nNull + 1
    }
  }
  transProbs = transProbs / rowSums(transProbs)
  return(transProbs)
}
computeEmissionProbs = function(rf) {
  # get the list of unque states
  states = sort(rf$classes)
  symbols = sort(rf$classes)
  S = length(states)
  # set up the transition Probability matrix
  emissionProbs = matrix(0, nrow=S, ncol=S)
  rownames(emissionProbs) = states
  colnames(emissionProbs) = states
  # loop through state sequence and count transitions
  for (k in 1:length(states)) {
    emissionProbs[k, ] = colSums(rf$votes[rf$groundTruth == states[k], ]) / sum(rf$groundTruth == states[k])
  }
  return(emissionProbs)
}
computePriorProbs = function(stateSeq) {
  states = sort(unique(stateSeq))
  states = states[!grepl("NULL", states)]
  priorProbs = tabulate(factor(stateSeq[stateSeq != "NULL"], levels=states))
  priorProbs = priorProbs / sum(priorProbs)
  return(priorProbs)
}