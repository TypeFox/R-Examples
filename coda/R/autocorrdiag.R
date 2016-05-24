### This is a replacement for the autocorrelation function
### which only gives the autocorrelations and not the cross correlations.

"autocorr.diag" <- function (mcmc.obj, ...) {
  UseMethod("autocorr.diag")
}

### Only looks at the diagonal elements of the autocorrlation matrix.
"autocorr.diag.mcmc" <- function (mcmc.obj,...) {
  ac <- autocorr(mcmc.obj,...)
  dd <- dim(ac)
  result <- matrix(NA,nrow=dd[1],ncol=dd[2],dimnames =dimnames(ac)[1:2])
  for (i in 1:dd[2]) {
    result[,i] <- ac[,i,i]
  }
  return (result)
}

### Looks at the average autocorrelation for all chains series.
"autocorr.diag.mcmc.list" <- function (mcmc.obj,...) {
  ac <- lapply(mcmc.obj,autocorr.diag.mcmc,...)
  result <- ac[[1]]
  if (length(ac) > 1) {
    for (chain in 2:length(ac)) {
      result <- result + ac[[chain]]
    }
  }
  return (result/length(ac))
}
