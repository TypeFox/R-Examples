estimateActivations <- function(cuesOutcomes, weightMatrix, unique=FALSE,...) {

  cues = rownames(weightMatrix)
  outcomes = colnames(weightMatrix)

  ## Check for NAs

  NA.cue_strings <- grep("(^NA_)|(_NA_)|(_NA$)",cuesOutcomes$Cues)
  NA.outcome_strings <- grep("(^NA)|(_NA_)|(_NA$)",cuesOutcomes$Outcomes)
  if(length(NA.cue_strings)>0)
    warning(paste("Potential NA's in ",length(NA.cue_strings)," 'Cues'.",sep=""))
  if(length(NA.outcome_strings)>0)
    warning(paste("Potential NA's in ",length(NA.outcome_strings)," 'Outcomes'.",sep=""))
  NA.cues <- which(is.na(cuesOutcomes$Cues))

  if("Outcomes" %in% names(cuesOutcomes))
    NA.outcomes <- which(is.na(cuesOutcomes$Outcomes))
  else
    {
      NA.outcomes = NULL
      warning("No 'Outcomes' column specified in 'cuesOutcomes'.")
    }

  if(length(NA.cues)>0)
    stop(paste("NA's in 'Cues': ",length(NA.cues)," cases.",sep=""))
  if(length(NA.outcomes)>0)
    stop(paste("NA's in 'Outcomes': ",length(NA.outcomes)," cases.",sep=""))

  obsCues = strsplit(as.character(cuesOutcomes$Cues), "_")
  uniqueObsCues = unique(unlist(obsCues))
  newCues = uniqueObsCues[!is.element(uniqueObsCues, cues)]

  if(length(newCues) > 0) {
    wnew = matrix(0, length(newCues), ncol(weightMatrix))
    rownames(wnew)=newCues
    colnames(wnew)=colnames(weightMatrix)
    w = rbind(weightMatrix, wnew)
    cues = c(cues, newCues)
  } else {
    w = weightMatrix
  }

  obsOutcomes = strsplit(as.character(cuesOutcomes$Outcomes), "_")
  uniqueObsOutcomes = unique(unlist(obsOutcomes))
  newOutcomes = uniqueObsOutcomes[!is.element(uniqueObsOutcomes, outcomes)]

  m = matrix(0, length(cues), nrow(cuesOutcomes))
  rownames(m) = cues
#  colnames(m) = cuesOutcomes$WordForm
#  colnames(m) = cuesOutcomes$Outcomes

  v = rep(0, length(cues))
  names(v) = cues

  for(i in 1:nrow(cuesOutcomes)) {
      v[obsCues[[i]]]=1
      m[,i] = v
      v[obsCues[[i]]]=0
  }

  a = t(w) %*% m

  if (unique) {
    activationMatrix <- unique(t(a))
  } else {
    activationMatrix <- t(a)
  }

  if (length(newCues)>0)
    warning(paste("There were ", length(newCues), " cues not present in 'weightMatrix'.",sep=""))  
  if (length(newOutcomes)>0)
    { # activationMatrix = cbind(activationMatrix,matrix(0,NROW(activationMatrix),length(newOutcomes),dimnames=list(NULL,newOutcomes)))
      warning(paste("There were ", length(newOutcomes), " outcomes not present in 'weightMatrix'.",sep=""))  
    }

  result <- list(activationMatrix = activationMatrix, newCues = newCues, newOutcomes = newOutcomes)
  
  return(result)

}
