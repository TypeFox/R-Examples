call.genes <-
function(CN, direction="gain", method.p.det="p_samples_amplicons", method.p.val="p_samples", n.validated=1, thres.p=0.05) {
  P <- attr(CN, "P")
  cn <- as.numeric(P[, "CN"]) 
  if (direction == "gain") ids.dir <- rownames(P)[which(cn > 2)]
  if (direction == "loss") ids.dir <- rownames(P)[which(cn < 2)] 
  amplicons <- rownames(CN) 
  genes <- sapply(strsplit(rownames(CN), "_|-"), function(x){x[1]})
  genes.u <- sort(unique(genes))
  ngenes <- length(genes.u)
  samples <- sort(colnames(CN))
  nsamples <- length(samples)
  calls <- matrix(nrow=ngenes, ncol=nsamples, "")
  rownames(calls) <- genes.u
  colnames(calls) <- samples
  CN.gene <- matrix(nrow=ngenes, ncol=nsamples)
  rownames(CN.gene) <- genes.u
  colnames(CN.gene) <- samples
  for (g in genes.u) {
    amplicons.g <- amplicons[which(genes == g)] 
    CN.gene[g, ] <- apply(CN[amplicons.g, , drop=FALSE], 2, mean)
    for (s in samples) {
      ids.gs <- paste(amplicons.g, s, sep="_")
      ids <- intersect(ids.gs, ids.dir)
      if (length(ids) > 0) {
        Q <- P[ids, , drop=FALSE]
        p <- as.numeric(Q[, "p"])
        p.det <- as.numeric(Q[, method.p.det])
        p.val <- as.numeric(Q[, method.p.val])
        index.min <- which.min(p)
        n.val <- length(which(p.val < thres.p)) - 1
        if (p.det[index.min] < thres.p && n.val >= n.validated) {
          mean.cn <- round(CN.gene[g, s], 1)
          min.p <- signif(p[index.min], 2)
          calls[g, s] <- paste(c("CN", "p", "nvalidated"), c(mean.cn, min.p, n.val), sep="=", collapse=", ")
        }
      }
    }
  }
  A <- matrix(nrow=ngenes, ncol=5, "")
  rownames(A) <- genes.u
  colnames(A) <- c("namplicons", "nsamples", "percent", "mean_CN", "samples")
  for (g in genes.u) {
    A[g, "namplicons"] <- length(which(genes == g))
    samples.g <- samples[which(calls[g, ] != "")]
    calls.g <- calls[g, samples.g]
    n <- length(samples.g)
    percent <- round(n/nsamples*100, 1)
    A[g, "nsamples"] <- n
    A[g, "percent"] <- percent
    if (length(samples.g) > 0) {
      cns <- CN.gene[g, samples.g]
      A[g, "mean_CN"] <- round(mean(cns), 1)
      if (direction == "gain") index.ord <- order(cns, decreasing=TRUE)
      if (direction == "loss") index.ord <- order(cns, decreasing=FALSE)
      A[g, "samples"] <- paste(samples.g[index.ord], " (", calls.g[index.ord], ")", sep="", collapse=", ")
    }
  }
  B <- matrix(nrow=nsamples, ncol=2)
  rownames(B) <- samples
  colnames(B) <- c("ngenes", "genes")
  for (s in samples) {
    genes.s <- genes.u[which(calls[, s] != "")]
    calls.s <- calls[genes.s, s]
    B[s, "ngenes"] <- length(genes.s)
    if (length(genes.s) > 0) {
      cns <- CN.gene[genes.s, s]
      if (direction == "gain") index.ord <- order(cns, decreasing=TRUE)
      if (direction == "loss") index.ord <- order(cns, decreasing=FALSE)
      B[s, "genes"] <- paste(genes.s[index.ord], " (", calls.s[index.ord], ")", sep="", collapse=", ")
    }
  }
  ncalls <- sum(as.numeric(A[, "nsamples"]))
  cat(ncalls, "/", ngenes*nsamples, " detected and validated genes with CN ", direction, " called using ", method.p.det, " for detection and ", method.p.val, " for validation\n", sep="")
  attr(calls, "CN") <- CN.gene
  attr(calls, "gene") <- A
  attr(calls, "sample") <- B
  return(calls)
}
