# library(RColorBrewer)

### Modified from q.map().
q.map.new <- function(x, filename, start = 1,
    levels = seq(0.01, 0.10, length.out = 10),
    margs = c(0.1,0.1,0.1,0.1), l.col = "red"){

  K <- dim(x)[1]
  lev.n <- length(levels) + 1

  colors <- RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")[9:1]
  colors <- c(colors[1:3], rgb(0.941, 0.204, 0.137),
              colors[4:6], rgb(0.996, 0.776, 0.380), colors[7:9])
  p <- 11
  col.mat <- x
  levels <- c(0, levels, 1)
  for (i in 1:p){
    ind <- (x >= levels[i]) & (x < levels[i+1])
    col.mat[ind] <- colors[i]
  }
  ind <- x >= levels[p]
  col.mat[ind] <- colors[i]

  par(mar = margs)
  plot( c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "", axes = FALSE)
  box()

  grid <- seq(0, 1, length.out = K + 2)

  label.tmp <- (1:(K + 1)) + start - 1
  label.text <- as.character(label.tmp)
#  label.text[label.tmp > 9] <-
#    c("a", "b", "c", "d", "e", "f", "g", "h", "i",
#      "j", "k", "l", "m", "n", "o", "p", "q", "r",
#      "s", "t", "u", "v", "w", "x", "y", "z")[label.tmp[label.tmp > 9] - 9]

  for (i in 1:K){
    for (j in 2:(K+1)){
      polygon(c(grid[j], grid[j+1], grid[j+1], grid[j]),
              c(grid[K-i+1], grid[K-i+1], grid[K-i+2], grid[K-i+2]),
              col = col.mat[i,j-1], border = NA)

      if (i == (j - 1)){
        x1 <- (grid[j-1] + grid[j]) / 2
        y1 <- (grid[K-i+1] + grid[K-i+2]) / 2
        text(x1, y1, label.text[i], cex = 1.5)
#        text(x1, y1, i+start-1, cex = 2)
        y1 <- (grid[K+1] + grid[K+2]) / 2
        x1 <- (grid[j] + grid[j+1]) / 2
        text(x1, y1, label.text[i + 1], cex = 1.5)
#        text(x1, y1, i+start, cex = 2)
      }

      lines(c(grid[i+1], 1), c(grid[K-i+1], grid[K-i+1]), col = l.col)
      lines(c(grid[K-i+2], grid[K-i+2]), c(grid[i], 1 - grid[2]), col = l.col)
    }
  }

  for (i in 1:K){
    for (j in 2:(K+1)){
      lines(c(grid[i+1], 1), c(grid[K-i+1], grid[K-i+1]), col = l.col)
      lines(c(grid[K-i+2], grid[K-i+2]), c(grid[i], 1 - grid[2]), col = l.col)
    }
  }

  lines(c(grid[2], 1), c(1 - grid[2], 1 - grid[2]), col = l.col)
  lines(c(1, 1), c(0, 1 - grid[2]), col = l.col)

  legend( 0, 0.5, legend = c("< 0.01", "0.01 - 0.02", "0.02 - 0.03",
                               "0.03 - 0.04", "0.04 - 0.05", "0.05 - 0.06",
                               "0.06 - 0.07", "0.07 - 0.08", "0.08 - 0.09",
                               "0.09 - 0.10", "> .10"),
         col = colors, pch = rep(15, p))
         # col = colors, pch = rep(15, p), cex = 1.3, box.lty = 0)
}


#### q-value
bh.fdr <- function(p){
  # This function computes q-values using Benjamini and Hochberg's (1995)
  # approach for controlling FDR.
  # Author: Dan Nettleton
  m <- length(p)
  k <- 1:m
  ord <- order(p)
  p[ord] <- (p[ord] * m) / (1:m)
  qval <- p
  for(i in (m - 1):1) {
    qval[ord[i]] <- min(c(qval[ord[i]], qval[ord[i + 1]]))
  }
  return(qval)
}

plotp <- function(pv.un, dim, start = 1, file.name = "un"){
  pv.m <- matrix(NA, ncol = dim, nrow = dim)
  pv.m[lower.tri(pv.m, diag = TRUE)] <- pv.un
  pv.m <- t(pv.m)
  q.map.new(pv.m, paste(file.name, "-p.pdf", sep = ""), start = start)
}
plotq <- function(pv.un, dim, start = 1, file.name = "un"){
  qv <- bh.fdr(pv.un)
  qv.m <- matrix(NA, ncol = dim, nrow = dim)
  qv.m[lower.tri(qv.m, diag = TRUE)] <- qv
  qv.m <- t(qv.m)
  q.map.new(qv.m, paste(file.name, "-p.pdf", sep = ""), start = start)
}
