IsoBHPlot <- function(rp, FDR, stat = c("E2", "Williams", "Marcus", "M", "ModifM")){

  stat <- match.arg(stat)
    
  Probe.ID <- rp[, 1]

  rpraw <- switch(stat,
      E2 = rp[, 2],
      Williams = rp[, 3],
      Marcus = rp[, 4],
      M = rp[, 5],
      ModifM = rp[, 6])

  #procs <- c("Bonferroni", "Holm", "BH", "BY")
  #res <- mt.rawp2adjp(rpraw, procs)  # from multtest
  #adjp <- res$adjp[order(res$index), ]

  adjp <- cbind(rpraw, p.adjust(rpraw, "BH"),p.adjust(rpraw,"BY"))
  

  plot(1:nrow(rp), sort(adjp[,1]),
       col = 4, pch = ".", lty = 1, xlab = "index",
       ylab = "Adjusted P values")

  lines(1:nrow(rp), sort(adjp[,1]), lty = 1, col = 1)
  lines(1:nrow(rp), sort(adjp[,2]), lty = 4, col = 2)
  lines(1:nrow(rp), sort(adjp[,3]), lty = 5, col = 3)
  abline(FDR, 0, lty = 6)

  legend(nrow(rp) / 2,
         0.3, col = c(1,2,3), c("Raw P","BH(FDR)","BY(FDR)"), lty = c(1,4:5))
  title(paste(stat, ": Adjusted p values by BH and BY", sep = ""))
}

