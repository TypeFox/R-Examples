print.lpcfdrobj <- function(x,...){
  lpcfdr.obj <- x
  cat("Call:\n",fill=F)
  dput(lpcfdr.obj$call)
  cat(fill=T)
  cat(c("Pi0:  ", round(lpcfdr.obj$pi0,4)),fill=T)
  cat(c("Soft Thresh:  ", round(lpcfdr.obj$soft.thresh,4)),fill=T)
  numGenesT <- NULL
  numGenesLPC <- NULL
  for(thr in c(.001,.01,.05,.1,.2,.4)){
    numGenesT <- c(numGenesT, sum(lpcfdr.obj$fdrt <= thr))
    numGenesLPC <- c(numGenesLPC, sum(lpcfdr.obj$fdrlpc <= thr))
  }
  cat(fill=T)
  cat("FDR Cut-off and # Genes with FDRs below that cut-off:", fill=T)
  mat <- cbind(numGenesT, numGenesLPC)
  dimnames(mat) <- list(c(.001,.01,.05,.1,.2,.4), c("# Genes with T FDR < Cutoff", "# Genes with LPC FDR < Cutoff"))
  print(mat, quote=FALSE)
  invisible()
}
