pamr.geneplot <- function(fit, data, threshold) {
  par(pch = 1, col = 1)
  geneid <- data$geneid
  if(is.null(geneid)) {
    geneid <- as.character(1:nrow(data$x))
  }
  if(is.null(fit$newy)) {
    y <- factor(data$y[fit$sample.subset])
  }
  else {
    y <- factor(fit$newy[fit$sample.subset])
  }
  x <- data$x[fit$gene.subset, fit$sample.subset]
  geneid <- geneid[fit$gene.subset]
  nc <- length(unique(y))
  aa <- pamr.predict(fit, x, threshold = threshold, type = "nonzero")
  cen <- pamr.predict(fit, x, threshold = threshold, type = "cen")
  d <- (cen - fit$centroid.overall)[aa,  ]/fit$sd[aa]
  oo <- order( - apply(abs(d), 1, max))
  aa <- aa[oo]
  ngenes <- length(aa)
  o <- order(y)
  xx <- x[aa, o]
  geneid <- geneid[aa]
  nc <- length(unique(y))
  nn <- c(0, cumsum(table(y)))
  nrow <- trunc(sqrt(ngenes)) + 1
  ncol <- trunc(sqrt(ngenes)) + 1
  if(nrow * (ncol - 1) >= ngenes) {
    ncol <- ncol - 1
  }
  par(mfrow = c(nrow, ncol))
  for(i in 1:ngenes) {
    plot(1:ncol(xx), xx[i,  ], type = "n", xlab = "sample", ylab = 
         "expression", axes = FALSE)
    box()
    axis(2)
    for(j in 1:nc) {
      j1 <- nn[j] + 1
      j2 <- nn[j] + table(y)[j]
      points(j1:j2, xx[i, j1:j2], col = j + 1)
    }
    title(main = as.character(geneid[i]))
    for(j in 1:(nc - 1)) {
      abline(v = cumsum(table(y))[j] + 0.5, lty = 2)
    }
    if(i == 1) {
      h <- c(0, table(y))
      for(j in 2:(nc + 1)) {
        text(sum(h[1:(j - 1)]) + 0.5 * h[j], max(xx[i,  
                                                    ]), label = levels(y)[j - 1], col = j)
      }
    }
  }
  par(mfrow = c(1, 1))
}


