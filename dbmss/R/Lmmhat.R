Lmmhat <-
function(X, r = NULL, ReferenceType = "", CheckArguments = TRUE) {
  
  if (CheckArguments)
    CheckdbmssArguments()
  
  KLmm <- Kmmhat(X, r, ReferenceType, CheckArguments)
  Columns <- names(KLmm)[-1]
  for(i in Columns) {
    KLmm[[i]] <- sqrt(KLmm[[i]]/pi)-KLmm$r
  }
  attr(KLmm, "ylab") <- attr(KLmm, "yexp") <- quote(L[mm](r))
  attr(KLmm, "fname") <- "L[mm]"
  return (KLmm)
}
