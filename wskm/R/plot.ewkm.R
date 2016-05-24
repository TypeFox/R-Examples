plot.ewkm <- function(x, ...)
{
  x  <- t(x$weights)
  rc <- rainbow(nrow(x), start=0, end=.3)
  cc <- rainbow(ncol(x), start=0, end=.3)
  # 120219 An alternative is to use pheatmap. Is there a ggplot2
  # alternative?
  hv <- heatmap(x, col = cm.colors(256), scale="column",
                RowSideColors = rc, ColSideColors = cc, margins=c(5,10),
                xlab = "Cluster", ...)
}
