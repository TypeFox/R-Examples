drulenew <- function(ruleset, modif, simutry, sig=0.95, data, noden,cutoff) 
{ 
  lik.chang <- rep(0, length(modif))
  vnam <- rep(NA, length(modif))
  colna <- colnames(ruleset)
  for (i in 1:length(modif))
  { 
    dvar <- order(modif, decreasing=T)[1]
    dff <- bound(noden)[dvar]
    if (modif[dvar] >= cutoff[dff])
    { 
      lik.chang[i] <- modif[dvar]
      if (ncol(ruleset)>1) vnam[i] <- names(modif)[dvar]
      else vnam[i] <- colna
      colna <- colnames(ruleset)[-dvar]   
      ruleset <- as.matrix(ruleset[,-dvar]) 
      colnames(ruleset) <- colna
      if (ncol(ruleset)>=1) modif <- change.log(ruleset, noden, simutry, data) - sum(lik.chang)
    }
    else break
  } 
  i <- sum(!is.na(vnam))
  if (i >= 1)
  {
    lik.change <- lik.chang[1:i]
    names(lik.change) <- vnam[1:i]
  }
  else lik.change <- -Inf
  return(lik.change)
}

