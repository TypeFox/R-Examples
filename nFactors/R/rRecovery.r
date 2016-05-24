"rRecovery" <-
function(R, loadings, diagCommunalities=FALSE) {
 recoveredR <- loadings %*% t(loadings)
 recovery   <- list(R = R, recoveredR = recoveredR, difference = R - recoveredR)
 if (diagCommunalities == FALSE) {diag(R)    <- NA; diag(recoveredR) <- NA }
 corr       <- cor(c(R),c(recoveredR), use="pairwise.complete.obs")
 recovery   <- list(recovery, cor = corr)
 return(recovery)
 }

 
