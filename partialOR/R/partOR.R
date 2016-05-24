# 12/31/2012 10:05:46 AM
partialOR <- function(dd,ci=0.95) {
  fit <- fitOR(dd)   
  reportOR(fit,dd,ci)
  return(invisible(NULL))
}

