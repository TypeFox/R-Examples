PMplot <- function(pmatrix, level=0.05, mtitle=NULL, xylabel=NULL, margin=5, legendx=0.73, newwd=TRUE){

  if (is.matrix(pmatrix)) {
    nr <- nrow(pmatrix)
    pmatrix[upper.tri(pmatrix, diag=TRUE)] <- NA
    if (is.null(rownames(pmatrix)))  rnpltm <- as.character(1:nrow(pmatrix))  else rnpltm <- rownames(pmatrix)
  }
  
  if (is.list(pmatrix)) {
    for (i in 1:length(pmatrix)) pmatrix[[i]][upper.tri(pmatrix[[i]], diag=TRUE)] <- NA
    pmatrix <- do.call(adiag, c(pmatrix, pad=NA))     
    nr <- nrow(pmatrix)
    if (is.null(xylabel)) rnpltm <- as.character(1:nrow(pmatrix)) else rnpltm <- xylabel
  }
  
  if (nr <= 3) {
    cat("\nThere is no plot for p-values matrix less than six values!\n")
  }else{
    if (is.null(mtitle)) mtitle <- paste("Level Plot of p-value Matrix")
    if (newwd) dev.new()
    pltmm <- t(pmatrix[nr:1,])
    if (level == 0.05) {
      pltm <- matrix(as.numeric(cut(as.numeric(pltmm), c(-0.1, 0.01, 0.05, 0.1, 1))), nrow=nr)
      pltmm <- matrix(as.numeric(droplevels(cut(as.numeric(pltmm), c(-0.1, 0.01, 0.05, 0.1, 1)))), nrow=nr)
    }else  pltm <- matrix(as.numeric(cut(as.numeric(pltmm), c(-0.1, level, 1))), nrow=nr)

    if (level == 0.05) pcolr <- c("#0D0DFF", "#5D5DFF", "#A1A1FF", "#E4E4FF") else pcolr <-  c("#0D0DFF", "#A1A1FF")
    colr <- pcolr[sort(unique(na.omit(as.numeric(pltm))))]
    max.len <- max(nchar(rnpltm))/6
    mar <- rep(margin, 2) #c(8, 8)
    op <- par(mar = c(mar[1] + max.len, mar[1] + max.len, 4, 4))
    zlim <- range(pltmm, na.rm=TRUE)
    image(pltmm, col = colr, axes = FALSE, main = mtitle, zlim = zlim)
    at1 <- (0:(nr - 1))/(nr - 1)
    tk <- at1 - 0.5/(nr - 1)
    if (max.len > 0.5) {
      axis(1, at = at1, labels = rnpltm, las = 2)
      axis(2, at = at1, labels = rnpltm[nr:1], las = 1)
    }else{
      axis(1, at = at1, labels = rnpltm)
      axis(2, at = at1, labels = rnpltm[nr:1], las = 1)
    }
    abline(h = tk[-1], v = tk[-1], col = "white")
    box()

    if (level == 0.05) {
      legen.lab <- c(expression(p > 0.1), expression("0.05 <  p " <=
        0.1), expression("0.01 < p " <= 0.05), expression(p <=
        0.01))[rev(5-sort(unique(na.omit(as.numeric(pltm)))))]
      legend(legendx, 0.99, legen.lab, pch = rep(15, length(colr)),
        col = rev(colr), pt.cex = 1.5, cex = 0.9)
    } else{
      legend(legendx, 0.99, title = paste("At", level,
      "level"), c("insignificant", "significant")[3-sort(unique(na.omit(as.numeric(pltm))))],
      pch = rep(15, length(colr)), col = rev(colr), pt.cex = 1.5, cex = 0.9)
    }
    par(op)
  }# end of if(nr <= 3)
}
