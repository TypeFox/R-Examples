levelplot.ewkm <- function(x, ...)
{
  x  <- x$weights
  dd.row <- as.dendrogram(hclust(dist(x)))
  row.ord <- order.dendrogram(dd.row)
 
  dd.col <- as.dendrogram(hclust(dist(t(x))))
  col.ord <- order.dendrogram(dd.col)

  my.colors <- function(n) rev(heat.colors(n))
  
  levelplot(x[row.ord, col.ord],
            aspect = "fill", pretty=TRUE,
            xlab="Cluster", ylab="",
            colorkey = list(space = "left", col=my.colors),
            col.regions=my.colors,
            legend =
            list(right =
                 list(fun = dendrogramGrob,
                      args =
                      list(x = dd.col, ord = col.ord,
                           side = "right",
                           size = 10)),
                 top =
                 list(fun = dendrogramGrob,
                      args =
                      list(x = dd.row, 
                           side = "top"))),
            ...)
}
