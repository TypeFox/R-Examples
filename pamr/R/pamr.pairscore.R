pamr.plotcen <- function(fit, data, threshold) {
  genenames <- data$genenames[fit$gene.subset]
  x <- data$x[fit$gene.subset, fit$sample.subset]
  clabs <- colnames(fit$centroids)
  scen <- pamr.predict(fit, data$x, threshold = threshold, type = "cent")
  dif <- scen - fit$centroid.overall
  nc <- length(unique(fit$y))
  o <- drop(abs(dif) %*% rep(1, nc)) > 0
  d <- dif[o,  ]
  nd <- sum(o)
  genenames <- genenames[o]
  xx <- x[o,  ]
  oo <- order(apply(abs(d), 1, max))
  d <- d[oo,  ]
  genenames <- genenames[oo]
  par(mar = c(1, 5, 1, 1), col = 1)
  plot(rep(2, nd) + d[, 1], 1:nd, xlim = c(0, 2*nc+1), ylim = c(1, nd + 3), 
       type = "n", xlab = "", ylab = "", axes = FALSE)
  box()
  abline(h = seq(nd), lty = 3, col = 7)
  jj <- rep(0, nd)
  for(j in 1:nc) {
    segments(jj + 2 * j, seq(nd), jj + 2 * j + d[, j], seq(nd), col
             = j + 1, lwd = 4)
    lines(c(2 * j, 2 * j), c(1, nd), col = j + 1)
    text(2 * j, nd + 2, label = clabs[j], col = j + 1)
  }
  g <- substring(genenames, 1, 20)
  text(rep(0, nd), seq(nd), label = g, cex = 0.4, adj = 0, col = 1)
}

