panel.bwplot.superpose <-
  function(x, y, ...,
           groups=groups,
           col=rep(trellis.par.get("superpose.symbol")$col, length=length(groups)),
           pch=trellis.par.get("box.dot")$pch,
           panel.groups=panel.bwplot.groups) {
    if (missing(groups)) stop("'groups=' argument must be specified.", call.=FALSE)
    par.old <- trellis.par.get()[c(
      "box.dot",
      "box.rectangle",
      "box.umbrella",
      "plot.symbol")]
    on.exit(trellis.par.set(par.old))
    panel=panel.superpose(as.position(x), y, ...,
      panel.groups=panel.groups,
      groups=groups,
      col=col,
      pch=pch)
  }

## panel.bwplot.groups <-
##     ## argument fill is captured and ignored
##   function(..., col, pch, fill, fill.alpha=NULL, group.number) {
##     if (is.null(col))
##       col <- group.number
##     else
##       col <- rep(col, length=group.number)[group.number]
##     trellis.par.set(
##       box.dot=list(col=col, pch=pch),
##       box.rectangle=list(col=col, alpha=ifelse(is.null(fill.alpha), 1, fill.alpha)),
##       box.umbrella=list(col=col, alpha=ifelse(is.null(fill.alpha), 1, fill.alpha)),
##       plot.symbol=list(col=col))
##     panel.bwplot(..., fill=ifelse(is.null(fill.alpha), 0, col))
##   }


panel.bwplot.groups <-
  function(..., col, pch, fill, fill.alpha=NULL, group.number) {
    ## argument fill is captured and ignored
    if (is.null(col))
      col <- group.number
    else
      col <- rep(col, length=group.number)[group.number]
    if (is.null(fill.alpha)) {
      trellis.par.set(
        box.dot=list(col=col, pch=pch),
        box.rectangle=list(col=col),
        box.umbrella=list(col=col),
        plot.symbol=list(col=col))
      panel.bwplot(..., fill=0)
    }
    else {
      trellis.par.set(
        box.dot=list(col=col, pch=pch),
        box.rectangle=list(col=col, alpha=fill.alpha),
        box.umbrella=list(col=col, alpha=fill.alpha),
        plot.symbol=list(col=col))
      panel.bwplot(..., fill=col)
    }
  }
