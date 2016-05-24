# output:
# list of significant genes
#    col 1: Probe.ID
#    col 2: row number
#    col 3: permutation p value
#    col 4: BH adjusted p value
#    col 5: q value
IsoTestSAM <- function(x, y, fudge, niter, FDR, stat) {
  qqstat <- Isoqqstat(x, y, fudge, niter)
  allfdr <- Isoallfdr(qqstat, , stat)
  del.table <- data.frame(allfdr)
  min_fdr <- min(na.exclude(del.table[, 5]))
  if (min_fdr > FDR) {
     FDR <- min_fdr
     delta <- min(na.exclude(del.table[del.table[,5] <= FDR, 1]))
     print("FDR cannot be obtained in this dataset")
  } else {
     delta <- min(na.exclude(del.table[del.table[,5] <= FDR, 1]))
  }
  qval <- Isoqval(delta, allfdr, qqstat, stat)
  q.value <- qval[[1]]

##adding the permutation p values and adjusted p values

  switch(stat,
      E2 = {
        qstat <- qqstat[[1]]
        dperm <- qqstat[[2]]},
      Williams = {
        qstat <- qqstat[[3]]
        dperm <- qqstat[[4]]},
      Marcus = {
        qstat <- qqstat[[5]]
        dperm <- qqstat[[6]]},
      M = {
        qstat <- qqstat[[7]]
        dperm <- qqstat[[8]]},
      ModifM = {
        qstat <- qqstat[[9]]
        dperm <- qqstat[[10]]
      })

  dperm.out <- sort(as.vector(dperm))
  d <- qstat[,1]

  p.value <- sapply(d,function(x) sum(abs(x) <=abs(dperm.out))/niter/length(d))
  
  #procs <- c("BH")
  #res <- mt.rawp2adjp(p.value, procs)
  #adj.p.value <- res$adjp[order(res$index), ]
  
  adj.p.value <- cbind(p.value, p.adjust(p.value, "BH"))
  

  qap <- cbind(q.value, adj.p.value)

  significant <- qap[qap[,3] <= FDR,, drop = FALSE] # TV: no drop for one row matrices
  sign.genes <- data.frame(row.names(y[significant[,1],]), significant)
  sign.genes1 <- sign.genes[order(significant[,2]),,drop=FALSE] # TV: no drop for one row matrices
  row.names(sign.genes1) <- 1:nrow(sign.genes1)
  names(sign.genes1) <- c("Probe.ID", "row.number","stat.val","qvalue","pvalue","adj.pvalue")

  return(list(sign.genes1, qqstat, allfdr))
}
