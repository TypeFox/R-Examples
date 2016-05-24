# plot with multi lines by group 

multilines <- function(XY, group=NULL, which=1:nf, sort=1, type='l', col=palette(), lwd=1, ...) {
  if (is.null(group)) group <- rep(1, nrow(XY))
  fact <- as.character(group)
  fact <- fact[!duplicated(fact)]
  nf <- length(fact)
  col <- rep(col, out.length=nf)
  lwd <- rep(lwd, out.length=nf)
  
  for (i in which) {
    xy <- subset(XY, subset=group==fact[i])
    if (sort %in% 1:2) {
      ord <- order(xy[, sort])
      xy <- xy[ ord, ]
    }
    lines(xy, type=type, col=col[i], lwd=lwd[i], ...)
  }
}
