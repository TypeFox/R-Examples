breslow <- function(f, riskset, entry, exit, trans, event){
  Ri <- riskset$Ri
  Ci <- riskset$Ci
  n <- length(entry)
  dummy <- rep(0, n)
  ef <- exp(f)
  for (j in 1:n){
    dummy[j] <- sum(ef[Ri[[j]]])
  }
  cbhr <- rep(0, n)
  dummy[which(dummy==0.0)] <- 1e-05
  for(i in 1:n){
    hi <- Ci[[i]]
    cbhr[i] <- sum(event[hi]/dummy[hi])
  }
  Q <- sort(unique(trans))
  A <- vector("list", length(Q))
  for(q in 1:length(Q)){
    hi <- which((trans == Q[q]) & (event == 1))
    A[[q]]$times <- exit[hi]
    A[[q]]$cbhr <- cbhr[hi]
    A[[q]]$cbhr <- A[[q]]$cbhr[order(A[[q]]$times)]
    A[[q]]$times <- A[[q]]$times[order(A[[q]]$times)]
  }
  names(A) <- paste("t", Q, sep="")
  return(A)}