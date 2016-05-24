qdc.zsc <- function(results) {
# Extract the scores of distinguishing statements from an object of class QmethodRes
  nstat <- results$brief$nstat
  nfactors <- results$brief$nfactors
  zsc <- results$zsc
  dac <- results$qdc$dist.and.cons
  names(dac) <- rownames(results$qdc)
  qdc.zsc <- matrix(NA, nrow=nrow(zsc), ncol=ncol(zsc), dimnames=dimnames(zsc))
  
  # Find which are the distinguishing statements
  for (i in 1:nfactors) {
    qdc.zsc[grep(paste0("all|f", i), dac), i] <- zsc[grep(paste0("all|f", i), dac), i]
  }
  return(qdc.zsc)
}
