runHMM <-
function(experiment, nstates=5, changeprob=0.0001, includeZeroState=T,errorRate=0.001){
  cat("Running HMM... \n")
  if(experiment$is_GC_normalized){readcounts <- experiment$CorrReadsprWindow}
  else{readcounts <- experiment$ReadsprWindow}
  newexp <- experiment
  cat("Analyzing chromosome", experiment$accession, "\n")
  chr_max <- max(readcounts)
  this_emission <- setupEmissionMatrix(nstates=nstates, mean=experiment$mean, variance=experiment$variance, absmax=chr_max, includeZeroState=includeZeroState,errorRate=errorRate)
  this_transition <- setupTransitionMatrix(nstates, remainprob=(1-changeprob), includeZeroState)
  cat("Finished setting up transition and emission matrices. Starting Viterbi algorithm...\n")
  copynumbers <- HMMcopynumber(readcounts, this_transition, this_emission, includeZeroState, experiment$windowlength, experiment$chrlength)
  cat("Finished running Viterbi algorithm. Assigning most probable states to individual segments...\n")
  # Now go through the genelist of the chromosome and insert Viterbi algorithm copy numbers.
  CN_HMM <- list()
    
  for (row in 1:nrow(experiment$genes)){
    start <- experiment$genes$Left[row]
    end <- experiment$genes$Right[row]
    CN_HMM[[row]] <- numeric()
    for (cnrow in 1:nrow(copynumbers)){
      state <- as.integer(copynumbers$State[cnrow])
      hmmstart <- as.integer(copynumbers$Startpos[cnrow])
      hmmend <- as.integer(copynumbers$Endpos[cnrow])
      within <- isbetween(c(start,end),hmmstart,hmmend)
      if(identical(within,c(T,T))){CN_HMM[[row]] <-state;next}
      if(identical(within,c(F,T)) | (identical(within, c(T,F)))) {CN_HMM[[row]] <- addState(CN_HMM[[row]],state);next}
      if (spans(c(start,end),hmmstart,hmmend)) {CN_HMM[[row]] <- addState(CN_HMM[[row]],state)}
    }
  }
  
  # Remove duplicate entries, NAs and sort states
  CN_HMM_RETURN <- vector()
  for(number in 1:length(CN_HMM)){
    toadd <- as.numeric(CN_HMM[number][[1]])
    toadd <- sort(toadd, na.last=NA)
    toadd <- unique(toadd)
    CN_HMM_RETURN[length(CN_HMM_RETURN)+1] <- paste(toadd,collapse=",")
  }
  newexp$genes$CN_HMM <- CN_HMM_RETURN
  newexp$HMMtable <- copynumbers
  cat("Copy number analysis complete. Use print to access results.\n")
  return(newexp)
}
