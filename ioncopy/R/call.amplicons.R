call.amplicons <-
function(CN, direction="gain", method.p="p_samples", thres.p=0.05) {
  amplicons <- sort(rownames(CN))
  namplicons <- length(amplicons)
  samples <- sort(colnames(CN))
  nsamples <- length(samples)
  P <- attr(CN, "P")
  p <- as.numeric(P[, method.p])
  cn <- as.numeric(P[, "CN"]) 
  index.p <- which(p < thres.p)
  if (direction == "gain") index.direction <- which(cn > 2)
  if (direction == "loss") index.direction <- which(cn < 2)
  index <- intersect(index.p, index.direction)
  ncalls <- length(index)
  cat(ncalls, "/", namplicons*nsamples, " amplicons with CN ", direction, " called using ", method.p, "\n", sep="")
  calls <- matrix(nrow=namplicons, ncol=nsamples, "")
  rownames(calls) <- amplicons
  colnames(calls) <- samples
  if (ncalls > 0) for (i in index) {
    cn.x <- round(cn[i], 1)
    p.x <- signif(as.numeric(P[i, "p"]), 2)
    a <- P[i, "amplicon"]
    s <- P[i, "sample"]
    calls[a, s] <- paste(c("CN", "p"), c(cn.x, p.x), sep="=", collapse=", ")
  }
  A <- matrix(nrow=namplicons, ncol=4, "")
  rownames(A) <- amplicons
  colnames(A) <- c("nsamples", "percent", "mean_CN", "samples")
  for (a in amplicons) {
    samples.a <- samples[which(calls[a, ] != "")]
    calls.a <- calls[a, samples.a]
    n <- length(samples.a)
    percent <- round(n/nsamples*100, 1)
    A[a, 1:2] <- c(n, percent)
    if (length(samples.a) > 0) {
      cns <- CN[a, samples.a]
      A[a, "mean_CN"] <- round(mean(cns), 1)
      if (direction == "gain") index.ord <- order(cns, decreasing=TRUE)
      if (direction == "loss") index.ord <- order(cns, decreasing=FALSE)
      A[a, "samples"] <- paste(samples.a[index.ord], " (", calls.a[index.ord], ")", sep="", collapse=", ")
    }
  }
  B <- matrix(nrow=nsamples, ncol=2, "")
  rownames(B) <- samples
  colnames(B) <- c("namplicons", "amplicons")
  for (s in samples) {
    amplicons.s <- amplicons[which(calls[, s] != "")]
    calls.s <- calls[amplicons.s, s]
    B[s, "namplicons"] <- length(amplicons.s)
    if (length(amplicons.s) > 0) {
      cns <- CN[amplicons.s, s]
      if (direction == "gain") index.ord <- order(cns, decreasing=TRUE)
      if (direction == "loss") index.ord <- order(cns, decreasing=FALSE)
      B[s, "amplicons"] <- paste(amplicons.s[index.ord], " (", calls.s[index.ord], ")", sep="", collapse=", ")
    }
  }
  attr(calls, "CN") <- CN
  attr(calls, "amplicon") <- A
  attr(calls, "sample") <- B
  return(calls)
}
