

rvquantile <- function(x, ...)
{
  UseMethod("rvquantile")
}

rvquantile.rv <- function(x, probs=c(0.025, 0.10, 0.25, 0.50, 0.75, 0.90, 0.975), ignoreInf=FALSE, ...)
{
  if (ignoreInf) {
    .f <- function (x) { quantile(x[is.finite(x)], probs=probs, ..., na.rm=TRUE) }
    t(rvsimapply(x, .f))
  } else {
    t(rvsimapply(x, quantile, probs=probs, ..., na.rm=TRUE))
  }
}

rvquantile.rvsummary <- function(x, probs=c(0.025, 0.10, 0.25, 0.50, 0.75, 0.90, 0.975), ...)
{
  Q <- t(sims(x))
  all_probs <- attr(Q, "quantiles")
  M <- NULL
  name <- character(0)
  # if (all(probs %in% all_probs)) ...
  for (p in probs) {
    ix <- (all_probs==p)
    if (any(ix)) {
      M <- cbind(M, Q[,ix,drop=FALSE])
    } else {
      name <- paste(p*100, "%", sep="")
      M <- cbind(M, NA)
      colnames(M)[ncol(M)] <- name
    }
  }
  return(M)
}

