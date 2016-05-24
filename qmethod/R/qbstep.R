qbstep <- function(subdata, subtarget, indet, nfactors, nqsorts, nstat, 
                   qmts=qmts, qmts_log=qmts_log, rotation="unknown", 
                   flagged=flagged, ...) {
  #--------------------------------------------------------------------
  # 1. Generate matrix of factor loadings
  cor.data <- cor(subdata, method="pearson")
  if (rotation=="unknown") rotation <- "none"
  loa <- as.data.frame(unclass(principal(cor.data, nfactors=nfactors, rotate=rotation, ...)$loadings))
  
  # Note (2015.12.17): the original line run 'principal()' directly: 
  # loa <- as.data.frame(unclass(principal(subdata, rotate="none", nfactors=nfactors)$loa))
  # However (funny enough!) principal() blocks the console when the data introduced are a square matrix (e.g. 30 observations and 30 statements); producing the correlation table first avoids that bug.
  #--------------------------------------------------------------------
  # 2. Apply solutions for indeterminacy issue of PCA bootstrap  
  if (indet == "none") {
    #loa <- as.data.frame(PCA(subdata, graph=FALSE)$var$coord[,c(1:nfactors)])
    loa <- as.matrix(unclass(varimax(as.matrix(loa))[[1]]))
  }
  if (indet == "procrustes") {
    #loa <- as.data.frame(PCA(subdata, graph=FALSE)$var$coord[,c(1:nfactors)])
    #caution: selecting rotation ="varimax" here implies that both varimax and Procrustes are used one on top of the other, and probably just one or the other should be used. For the qindtest though, the selected rotation is used
    procrustes <- qpcrustes(loa=loa, target=subtarget, nfactors=nfactors)
    loa <- procrustes
  }
  if (indet == "qindtest" | indet == "both") {
    loa <- as.matrix(unclass(varimax(as.matrix(loa))[[1]]))
    qindeterminacy <- qindtest(loa=loa, target=subtarget, 
                               nfactors=nfactors)
    loa <- qindeterminacy[[1]]
    if (indet == "both") {
      loa <- qpcrustes(loa=loa, target=subtarget, nfactors=nfactors)
    }
  }
  #--------------------------------------------------------------------
  # 3. Calculate z-scores and factor scores with the indeterminacy corrected factor loadings 'loa'
  flagged <- qflag(nstat=nstat, loa=loa)
  qstep <- qzscores(subdata, nfactors=nfactors,  
                    flagged=flagged, loa=loa, ...)
  #--------------------------------------------------------------------
  # 4. Export necessary results
  step_res <- list()
  step_res$flagged  <- list()
  step_res$zsc <- list()
  step_res$loadings <- list()
  
  qstep$flagged <- as.data.frame(qstep$flagged)
  qstep$loa     <- as.data.frame(qstep$loa)
  
  n <- 1
  while (n <= nfactors) {
    # Flagged q sorts
    step_res$flagged[n]  <- qstep$flagged[n]   #to append in qmbr[[n]][[1]]
    # z-scores
    step_res$zsc[n]      <- qstep$zsc[n]       #to append in qmbr[[n]][[2]]
    # Factor loadings
    step_res$loadings[n] <- qstep$loa[n]       #to append in qmbr[[n]][[3]]
    n <- n + 1
  }
  if (indet == "qindtest" | indet == "both") {
    qindt_log <- qindeterminacy[[2]]
    qindt     <- qindeterminacy[[3]]
    # Test results (logical)
    step_res$torder_res_log <- qindt[1]        #to append in qmts[1]
    step_res$tsign_res_log  <- qindt[2]        #to append in qmts[2]
    # Reports of solution implementation
    step_res$torder_res_report <- qindt_log[1] #to append in qmts_log[1]
    step_res$tsign_res_report  <- qindt_log[2] #to append in qmts_log[2]
  }
  return(step_res)
}