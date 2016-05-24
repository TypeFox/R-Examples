IsoSAMPlot <- function(qqstat, allfdr, FDR, stat){
  
  FDR50 <- allfdr[,5]
  FDR90 <- allfdr[,6]
  delta <- allfdr[,1]
  fp50 <- allfdr[,2]
  fp90 <- allfdr[,3]
  signnum <- allfdr[,4]

  par(mfrow = c(2, 2))

  # FDR vs. Delta
  plot(delta, FDR50, pch = ".", ylab = "FDR")
  lines(delta, FDR90, lty = 1)
  lines(delta, FDR50, lty = 2)
  abline(0.05, 0)
  abline(0.1, 0)
  legend(max(delta)*(1/2), 0.6, c("FDR90%","FDR50%"), lty = 1:2)
  title("a: plot of FDR vs. Delta")

  # sign vs. Delta
  plot(delta, signnum, pch = ".", ylab = "# of significant genes")
  lines(delta,signnum)
  title("b: plot of # of sign genes vs. Delta")

  # Fp vs. Delta
  plot(delta, fp50, pch = ".", ylab = "# of false positives")
  lines(delta, fp90, lty = 1)
  lines(delta, fp50, lty = 2)
  legend(max(delta)*(1/2), nrow(qqstat[[1]])/2, c("FP90%", "FP50%"), lty = 1:2)
  title("c: plot of # of FP vs. Delta")

  # obs vs. exp
  observed <- switch(stat,
      E2 = qqstat[[1]][, 1],
      Williams = qqstat[[3]][, 1],
      Marcus = qqstat[[5]][, 1],
      M = qqstat[[7]][, 1],
      ModifM = qqstat[[9]][, 1])
  
  expected <- switch(stat,
      E2 = qqstat[[1]][, 2],
      Williams = qqstat[[3]][, 2],
      Marcus = qqstat[[5]][, 2],
      M = qqstat[[7]][, 2],
      ModifM = qqstat[[9]][, 2])
  
  del.table <- data.frame(allfdr)
  min_fdr <- min(na.exclude(del.table[, 5]))
  if (min_fdr > FDR) {
    FDR <- min_fdr
    delta <- min(na.exclude(del.table[del.table[,5] <= FDR,1]))
    print("FDR cannot be obtained in this dataset")
  } else {
    delta <- min(na.exclude(del.table[del.table[,5] <= FDR,1]))
  }
                                            
  plot(expected, observed)
  abline(0, 1)
  abline(delta, 1, lty = 5, col = "blue")
  abline(-delta, 1, lty = 5, col = "blue")
   
  oe <- order(expected)
  q.mat <- switch(stat,
      E2 = qqstat[[1]][oe, ],
      Williams = qqstat[[3]][oe, ],
      Marcus = qqstat[[5]][oe, ],
      M = qqstat[[7]][oe, ],
      ModifM = qqstat[[9]][oe, ])
   
  if (sum(q.mat[,3] >= delta)>0) {
    x.exp <- min(q.mat[q.mat[,3] >= delta,2])
    # y.obs <- q.mat[which(q.mat[,2] == x.exp), 1]
    points(q.mat[q.mat[,2] >= x.exp,2], q.mat[q.mat[,2] >= x.exp,1], col = "red") 
  }

  if (sum(q.mat[,3] <= -delta)>0) {
    x.exp <- max(q.mat[q.mat[,3] <= -delta ,2])
    # y.obs <- q.mat[which(q.mat[,2] == x.exp), 1]
    points(q.mat[q.mat[,2] <= x.exp,2], q.mat[q.mat[,2] <= x.exp,1], col = "red") 
  }
  
  legend(min(expected), max(observed), paste("delta=", delta, sep=""), lty = 5)
  title("d: observed vs. expected statistics")
}
