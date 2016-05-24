print.tfdrobj <- function(x,...){
  tfdr.obj <- x
  cat("Call:\n",fill=F)
  dput(tfdr.obj$call)
  cat(fill=T)
  cat(c("Pi0:  ", round(tfdr.obj$pi0,4)),fill=T)
  numGenesT <- NULL
  for(thr in c(.001,.01,.05,.1,.2,.4)){
    numGenesT <- c(numGenesT, sum(tfdr.obj$fdrt <= thr))
  }
  cat(fill=T)
  cat("FDR Cut-off and # Genes with FDRs below that cut-off:", fill=T)
  mat <- matrix(numGenesT, ncol=1)
  dimnames(mat) <- list(c(.001,.01,.05,.1,.2,.4), c("# Genes with T FDR < Cutoff"))
  print(mat, quote=FALSE)
  invisible()
}
