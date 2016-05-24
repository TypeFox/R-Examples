wquantile <- function(wt=rep(1,length(x)), x, probs, already.sorted=FALSE, already.normalized=FALSE) {
  if(!already.sorted) {
    wt <- wt[o<-order(x)]
    x <- x[o]
  }
  if(!already.normalized) {
    wt <- wt/sum(wt)
  }
  x[findInterval(probs,cumsum(wt))]
}

wIQR <- function(wt=rep(1,length(x)), x, already.sorted=FALSE, already.normalized=FALSE) {
  diff(wquantile(wt, x, c(.25,.75), already.sorted, already.normalized))
}



