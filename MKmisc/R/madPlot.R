## Modifikation of function plot.cor of package "sma"
madPlot <- function (x, new = FALSE, col, maxMAD = 3, labels = FALSE, 
                     labcols = "black", title = "", protocol = FALSE, ...){
    n <- ncol(x)
    MAD <- x

    if(new) MAD <- madMatrix(x)

    if(missing(col))
      col <- rev(colorRampPalette(brewer.pal(8, "RdYlGn"))(128))

    if(protocol)
      layout(matrix(c(1, 2), 1, 2), widths = c(10, 3))
    else
      layout(matrix(c(1, 2), 1, 2), widths = c(10, 2.5))
    if(max(MAD) < maxMAD) 
      col.nr <- trunc(round(max(MAD), 2)/maxMAD*length(col))
    else 
      col.nr <- 1
    image(1:n, 1:n, MAD[, n:1], col = col[1:col.nr], axes = FALSE, 
          xlab = "sample index", ylab = "", ...)

    if (length(labcols) == 1) {
        axis(2, at = n:1, labels = labels, las = 2, cex.axis = 0.8,
            col.axis = labcols)
        axis(1, at = 1:n, labels = 1:n, las = 1, cex.axis = 0.8,
            col.axis = labcols)
    }

    if (length(labcols) == n) {
        cols <- unique(labcols)
        for (i in 1:length(cols)) {
            which <- (1:n)[labcols == cols[i]]
            axis(2, at = (n:1)[which], labels = labels[which],
                las = 2, cex.axis = 0.8, col.axis = cols[i])
            axis(1, at = which, labels = labels[which], las = 2,
                cex.axis = 0.8, col.axis = cols[i])
        }
    }
    title(title)
    box()

    x.bar <- seq(0, max(maxMAD, max(MAD, na.rm = TRUE)), length = length(col))
    x.small <- seq(x.bar[1], x.bar[length(x.bar)], length = 10)
    par(mar = c(5.1, 1, 4.1, 5))
    if(protocol){
      image(1, x.bar, matrix(x.bar, 1, length(x.bar)), axes = FALSE, xlab = "", ylab = "", col = col, ...)
      box()
      axis(4, at = c(x.small[2], x.small[9]), labels = c("dissimilar", "similar "), las = 2)
    }else{
      image(1, x.bar, matrix(x.bar, 1, length(x.bar)), axes = FALSE, xlab = "", ylab = "", col = col, ...)
      box()
      x.small <- seq(x.bar[1], x.bar[length(x.bar)], length = 10)
      Labels <- signif(rev(x.small), 2)      
      axis(4, at = rev(x.small), labels = Labels, las = 1)
    }
    
    layout(1)
    par(mar = c(5, 4, 4, 2) + 0.1)

    invisible()
}

