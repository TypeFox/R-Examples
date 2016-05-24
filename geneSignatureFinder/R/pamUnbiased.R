pamUnbiased <-
function(ddata) {
  n <- length(ddata)
  missing <- is.na(ddata)
  
  wdTmp <- ddata[!missing]
  nn <- n - sum(missing)
  
  ansClusters <- rep(NA, n)
  i <- 1
  while(i <= nn) {
    
    nMissing <- sum(missing)
    tmp <- pam(wdTmp[-i], 2)
    
    if(tmp$medoids[1] == tmp$medoids[2]) {
      ansClusters <- rep(0, n)
      missing <- rep(TRUE, n)
      break
    }
    
    clusters <- tmp$clustering
    # qui viene aggiornato l'elengo dei valori mancanti
    # con le osservazioni che fanno fallire PAM (valori anomali che generano cluster con 1 osservazione)
    ckTmp <- rep(FALSE, n)
    if(sum(ckTmp[!missing][-i] <- clusters == 1) < 0.025 * n)
      missing[!missing][-i] <- (missing[!missing][-i] + ckTmp[!missing][-i]) > 0
    
    ckTmp <- rep(FALSE, n)
    if(sum(ckTmp[!missing][-i] <- clusters == 2) < 0.025 * n)
      missing[!missing][-i] <- (missing[!missing][-i] + ckTmp[!missing][-i]) > 0
    
    if(nMissing < sum(missing)) {
      nMissing <- sum(missing)
      wdTmp <- ddata[!missing]
      nn <- n - sum(missing)
      ansClusters <- rep(NA, n)
      i <- 1
      next
    }
    
    m1 <- median(wdTmp[-i][clusters == 1])
    m2 <- median(wdTmp[-i][clusters == 2])
    
    s1 <- bw.nrd0(wdTmp[-i][clusters == 1])
    s2 <- bw.nrd0(wdTmp[-i][clusters == 2])
    # classificazione dell'osservazione in base alla regola di min distanza (standardizzata)
    ansClusters[!missing][i] <- which.min(abs((wdTmp[i] - c(m1, m2))/c(s1, s2)))
    i <- i + 1
  }
    
  tmp <- sum(ansClusters[!missing] == 1)
  if(min(tmp, nn - tmp)/nn < 0.025) {
    ansClusters <- rep(0, n)
    missing <- rep(TRUE, n)
  }
  ans <- list()
  ans$clusters <- ansClusters
  ans$missing  <- missing
  return(ans)
}
