plot.cover <-
function(x, shape = c("square", "circle", "equispaced"), pch = 21, cex = max(nchar(rownames(x))) + 2, bg = "white", ...) {
  V <- -vertices(x, shape = shape[1])
  xlim <- c(min(V$x), max(V$x))*1.3
  ylim <- c(min(V$y), max(V$y))*1.3
  plot(V, panel.first=drawedges(x, V), axes=FALSE, xlab="", ylab="", xlim=xlim, ylim=ylim, pch = pch, cex = cex, bg = bg, ...)
  text(V, labels = rownames(x), cex = 0.75)
}
