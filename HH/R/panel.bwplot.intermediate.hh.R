"panel.bwplot.intermediate.hh" <-
function (x, y,
          horizontal = TRUE,
          transpose=!horizontal,
          pch,
          col,
          at,  ## formerly S-Plus only, now totally ignored
          ...
          )
{
  if (missing(horizontal) && !missing(transpose))
    horizontal <- !transpose

  fac.levels <- if (horizontal) levels(y) else levels(x)
  box.par <- list(box.dot=trellis.par.get("box.dot"),
                  box.rectangle=trellis.par.get("box.rectangle"),
                  box.umbrella=trellis.par.get("box.umbrella"),
                  plot.symbol=trellis.par.get("plot.symbol"))
  old.box.par <- box.par
  on.exit(trellis.par.set(old.box.par))
  tpg <- trellis.par.get("superpose.line")
  tpg.col <- rep(tpg$col, length=length(fac.levels))
  if (!missing(pch)) pch <- rep(pch, length=length(fac.levels))
  if (!missing(col)) tpg.col <- rep(col, length=length(fac.levels))

  for (i in seq(along=fac.levels)) {
    if (!missing(pch)) {
      box.par$box.dot$pch <- pch[i]
      box.par$plot.symbol$pch <- pch[i]
    }
    for (j in names(box.par)) {
      box.par[[j]]$col <- tpg.col[i]
      trellis.par.set(j, box.par[[j]])
    }

    if (horizontal) {
      ii <- as.position(y[y == fac.levels[i]])
      xy <- x[y == fac.levels[i]]
      panel.bwplot(xy, ii, horizontal=horizontal, ...)
    }
    else {
      yx <- y[x == fac.levels[i]]
      ii <- as.position(x[x == fac.levels[i]])
      panel.bwplot(ii, yx, horizontal=horizontal, ...)
    }
  }
}
