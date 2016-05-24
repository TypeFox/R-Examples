print.lpcobj <- function(x,...){
  lpc.obj <- x
  cat("Call:\n",fill=F)
  dput(lpc.obj$call)
  cat(fill=T)
  cat(c("Soft Thresh:  ", round(lpc.obj$soft.thresh,4)),fill=T)
  cat(c("Number of non-zero coefficients:   ", sum(lpc.obj$coefs!=0)), fill=T)
  cat(fill=T)
  if(sum(lpc.obj$coefs!=0)>0){
    cat("Information about non-zero coefficients:", fill=T)
    mat <- cbind(which(lpc.obj$coefs!=0), round(lpc.obj$coefs[lpc.obj$coefs!=0],4))
    dimnames(mat) <- list(1:sum(lpc.obj$coefs!=0), c("Coefficient Index", "Coefficient Value"))
    print(mat, quote=FALSE)
  }
  invisible()
}
