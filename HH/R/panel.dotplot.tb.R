panel.dotplot.tb <- function(x, y, factor=.1,
                             jitter.data=TRUE, horizontal=TRUE,
                             max.freq=max(sapply(subsets, length)),
                             ...) {
  
  ## current.panel.limits()

  o.y.x <- order(y, x)
  x <- x[o.y.x]
  y <- y[o.y.x]

  present <- !(is.na(x) | is.na(y))
  x <- x[present]
  y <- y[present]
  
  subsets <- tapply(y, list(y,x), c)

  incr <- factor/max.freq
  
  y.j <- unlist(sapply(t(subsets),
                       function(x, incr) x+(seq(along=x)-1)*incr,
                       incr=incr))
  y.j <- y.j[!is.na(y.j)]


  if.R(r={
    panel.dotplot(x, y.j, levels.fos=unique(y), ...)
  },s={
    dot.symbol <- trellis.par.get("dot.symbol")
    dot.line <- trellis.par.get("dot.line")
    panel.abline(h = unique(y),
                 lwd = dot.line$lwd,
                 lty = dot.line$lty,
                 col = dot.line$col)
    ## In S-Plus 8, pch col cex font in the list(...) takes precedence
    ## because it appears second.  There is no complaint over multiple
    ## uses of these arguments.
    points(x, y.j,
           pch = dot.symbol$pch, col = dot.symbol$col,
           cex = dot.symbol$cex, font = dot.symbol$font,
           ...)
  })
}
