"panel.interaction2wt" <-
function(x, y, subscripts,
         responselab, trace.values,
         factor.levels, factor.position,
         fun=mean,
         se,
         type="l",
         ...,
         box.ratio,
         simple=FALSE,
         simple.offset,
         simple.scale,
         simple.pch,
         data.x,
         col.by.row=TRUE,
         key.in=NULL ## list of key arguments for S-Plus
         ) {

  if ("o" %in% type || "b" %in% type)
    type <- c(type, "p", "l")


  if (simple) {
    pch.name <- names(data.x)[3-current.row()]
    other.name <- names(data.x)[current.row()]
    pch <- list(NULL, NULL)
    names(pch) <- c(pch.name, other.name)
    lengths <- sapply(factor.levels, length)
    pch[[other.name]] <- rep(simple.pch[[pch.name]], lengths[other.name])
    pch[[pch.name]] <- rep(simple.pch[[other.name]], each=lengths[pch.name])
  }
  ## if.R(r={
    tcL <- trellis.currentLayout()
    rows.per.page <- dim(tcL)[1]
    cols.per.page <- dim(tcL)[2]
    cell <- panel.number()
    row.panel <- row(tcL)[tcL==cell]
    column.panel <- col(tcL)[tcL==cell]

    these.labels <- c(trace.factor=names(factor.levels)[row.panel],
                      x.factor=names(factor.levels)[column.panel])
  ## },s={
  ##   rows.per.page <- length(factor.levels)
  ##   cols.per.page <- length(factor.levels)
  ##   cell <- get("cell", frame=sys.parent())

  ##   tcL <- matrix(seq(length=rows.per.page*cols.per.page),
  ##                 rows.per.page, cols.per.page,
  ##                 byrow=TRUE)

  ##   row.panel <- row(tcL)[tcL==cell]
  ##   column.panel <- col(tcL)[tcL==cell]

  ##   these.labels <- c(trace.factor=names(factor.levels)[row.panel],
  ##                     x.factor=names(factor.levels)[column.panel])
  ## }
  ##      )

  trace.name <- these.labels["trace.factor"]
  trace.position <- factor.position[[trace.name]]
  trace.levels <- factor.levels[[trace.name]]
  x.name <- these.labels["x.factor"]
  x.position <- factor.position[[x.name]]
  x.levels <- factor.levels[[x.name]]

  tpg <- trellis.par.get("superpose.line")
  tpg.col <- rep(tpg$col, length=length(trace.levels))
  tpg.lty <- rep(tpg$lty, length=length(trace.levels))
  tpg.lwd <- rep(tpg$lwd, length=length(trace.levels))
  tpg <- trellis.par.get("superpose.symbol")
  tpg.cex <- rep(tpg$cex, length=length(trace.levels))
  tpg.pch <- rep(tpg$pch, length=length(trace.levels))
  tpg.alpha <- rep(tpg$alpha, length=length(trace.levels))
  tpg.fill <- rep(tpg$fill, length=length(trace.levels))
  tpg.font <- rep(tpg$font, length=length(trace.levels))

  ## panels
  if (trace.name ==  x.name) {## main diagonal
    if (simple && cols.per.page==2) { ## simple effects
      other.name <- names(factor.levels)[names(factor.levels) != x.name]
      other.factor <- data.x[, other.name]
      x.factor <- data.x[, x.name]
      ioh.list <- list(x.factor, other.factor)
      names(ioh.list) <- c(x.name, other.name)
      if (!missing(simple.offset))
        ioh.list$b.offset <- simple.offset[[other.name]]
      if (!missing(simple.scale))
        ioh.list$b.scale=simple.scale[[other.name]]
      x.simple <- do.call("interaction.positioned", ioh.list)
## recover()
      tpg.col.simple <-
        if (col.by.row)
          tpg.col
        else
          rep(tpg$col, length=length(levels(other.factor)))

      col.subscripts <-
        if (col.by.row)
          rep(seq(length(x.levels)), each=length(levels(other.factor)))
        else
          rep(seq(length(levels(other.factor))), length(x.levels))

      if.R(r=
           panel.bwplot.intermediate.hh(x.simple, y,
                                        horizontal=FALSE,
                                        col=tpg.col.simple[col.subscripts],
                                        pch=pch[[x.name]],
                                        box.ratio=box.ratio,
                                        ...)
           ,s=
           panel.bwplot.intermediate.hh(as.numeric(x.simple), y,
                                        horizontal=FALSE,
                                        col=tpg.col.simple[col.subscripts],
                                        pch=pch[[x.name]],
                                        box.ratio=box.ratio,
                                        ...)
           )
    }
    else { ## marginal main effects
      x.factor <- data.x[, x.name]
      if.R(r=
           panel.bwplot.intermediate.hh(x.factor, y,
                                        horizontal=FALSE,
                                        box.ratio=box.ratio,
                                        ...)
           ,s=
           panel.bwplot.intermediate.hh(as.numeric(x.factor), y,
                                        horizontal=FALSE,
                                        box.ratio=box.ratio,
                                        ...)
           )
    }
  }
  else { ## off-diagonal
    if (simple) {
      suff.data <-
        sufficient(data.frame(y, x, tt=trace.values[subscripts]),
                   yname="y", c("x", "tt"))

      ioh.list <- list(suff.data$x, suff.data$tt)
      names(ioh.list) <- c(x.name, trace.name)
      position(ioh.list[[x.name]]) <- x.position
      position(ioh.list[[trace.name]]) <- trace.position
      if (!missing(simple.offset))
        ioh.list$b.offset <- simple.offset[[trace.name]]
      if (!missing(simple.scale))
        ioh.list$b.scale=simple.scale[[trace.name]]
      x.simple <- do.call("interaction.positioned", ioh.list)
##browser()
      if (missing(se) || (is.logical(se) && !se))
        panel.intxplot(y=suff.data$y,
                       x=x.simple,
                       subscripts=seq(length(suff.data$y)),
                       groups=suff.data$tt,
                       offset.use=FALSE,
                       rug.use=TRUE,
                       type=type,
                       ...)
      else {
        if (is.logical(se))
          suff.data.se <- suff.data$sd/sqrt(suff.data$nobs)
        else {
          suff.data.se <- eval(se, suff.data)
        }
        panel.intxplot(y=suff.data$y,
                       x=x.simple,
                       subscripts=seq(length(suff.data$y)),
                       groups=suff.data$tt,
                       offset.use=FALSE,
                       offset=position(x.simple), ##simple.offset[[trace.name]],
                       rug.use=TRUE,
                       se=suff.data.se,
                       type=type,
                       ...)
      }
    }
    else {
      tab <- tapply(y, list(x, trace.values[subscripts]), fun)
      if.R(r={}, s=panel.lines <- lines)
      for (j in 1:ncol(tab)) {
        panel.lines(x=x.position, y=tab[,j], col=tpg.col[j], lty=tpg.lty[j], lwd=tpg.lwd[j])
        if ("p" %in% type)
          panel.points(x=x.position, y=tab[,j], col=tpg.col[j], pch=tpg.pch[j], cex=tpg.cex[j])
      }
    }
  }



  ## x labels
  if (row.panel==1) {
    if.R(r={},
         s={
           axis(1, at=x.position, labels=x.levels)
           mtext(x.name, side=1, line=2, at=mean(par()$usr[1:2]))
          })
}

  ## trace key
  if (column.panel==1) {

    key.list <- list(title=trace.name,  ## S-Plus only
                     cex.title=1,
                     corner=c(0,.5), border=TRUE,
                     text=list(text=trace.levels, cex=.8),
                     lines=list(col=tpg.col, lty=tpg.lty, lwd=tpg.lwd))
    key.list[names(key.in)] <- key.in
    if.R(r={},
         s={
           if (is.null(key.list$x))
             key.list$x <- par()$usr[1]-.6*diff(par()$usr[1:2])
           do.call("key", key.list)
         }
         )
  }


  ## y label
  if (column.panel==cols.per.page) {
    if.R(r={},
         s={
           mtext(responselab, side=4, line=3, at=mean(par()$usr[3:4]))
         }
         )
  }
}


axis.i2wt <-
  function(side, scales, ...)
{
  axis.options <- lattice.getOption("axis.options")

  tcL <- trellis.currentLayout()
  rows.per.page <- dim(tcL)[1]
  cols.per.page <- dim(tcL)[2]
  cell <- panel.number()
  row.panel <- row(tcL)[tcL==cell]
  column.panel <- col(tcL)[tcL==cell]

  if (side != "right")
  axis.default(side = side, scales = scales, ...)

  if (side == "bottom")
    {
      if (row.panel == 1)
        {

          axis.options <- axis.options$bottom
          at2     <- axis.options$at2
          labels2 <- axis.options$labels2
          rot2    <- axis.options$rot2
          ticks2  <- axis.options$ticks2
          at3     <- axis.options$at3
          labels3 <- axis.options$labels3
          rot3    <- axis.options$rot3
          ticks3  <- axis.options$ticks3
          if (is.null(rot2)) rot2 <- 0
          if (is.null(rot3)) rot3 <- 0
          if (is.null(ticks2)) ticks2 <- !is.null(labels2)
          if (is.null(ticks3)) ticks3 <- FALSE

          axis.units <- lattice.getOption("axis.units")$outer$bottom

          ## space occupied by standard set of ticks
          tick.units <-
            ## if (is.null(scales$at)) unit(0, "lines")
            ## else
              with(axis.units$tick, unit(scales$tck * x, units))

          ## space occupied by standard set of labels, assume 1.2
          ## lines, but we have enough information to compute
          ## explicitly (see how 'lab.grob' is computed in layout.R)

          label.units <-
            if (is.null(scales$at) || !scales$at) unit(0, "lines")
            else unit(1.2, "lines")

          ## want second set of ticks to extend beyond this.

          if (!is.null(at2))
            at2i <-
              if (is.list(at2))
                at2[[column.panel]]
              else at2

          if (!is.null(labels2))
            labels2i <-
              if (is.list(labels2))
                labels2[[column.panel]]
              else labels2

          labels2.units <-
            if (!is.null(at2) && !is.null(labels2))
              unit(1, "grobheight",
                   data = list(textGrob(unlist(labels2),
                     rot = rot2,
                     x=seq(along=unlist(labels2)),
                     just = if (rot2==90)
                     c("right", "center")
                     else c("center", "bottom"),
                     gp=do.call("gpar", trellis.par.get()$axis.text))))
            else
              unit(0, "lines")
          labels2.units <- convertHeight(labels2.units, "mm")

          if (ticks2 && !is.null(at2))
            grid.segments(x0 = unit(at2i, "native"),
                          x1 = unit(at2i, "native"),
                          y0 = unit(0, "npc"),
                          y1 = unit(0, "npc") - tick.units -
                          label.units)

          if (!is.null(at2) && !is.null(labels2))
            grid.text(labels2i, rot = rot2,
                      x = unit(at2i, "native"),
                      y = unit(0, "npc") - tick.units - label.units -
                      unit(1, "mm"),
                      just = if (rot2==90) c("right", "center") else
                      c("center", "top"),
                      gp=do.call("gpar", trellis.par.get()$axis.text))

          ## other set of labels (without ticks)

          if (ticks3 && !is.null(at3))
            grid.segments(x0 = unit(at3, "native"),
                          x1 = unit(at3, "native"),
                          y0 = unit(0, "npc"),
                          y1 = unit(0, "npc") - tick.units -
                          label.units - unit(1, "mm"))

          if (!is.null(at3) && !is.null(labels3))
            grid.text(labels3, rot=rot3,
                      x = unit(at3, "native"),
                      y = unit(0, "npc") - tick.units - label.units -
                      labels2.units - unit(8, "mm"),
                      just = c("center", "center"),
                      gp=do.call("gpar", trellis.par.get()$par.xlab.text))
          if (is.null(at3) && !is.null(labels3))
            grid.text(labels3[column.panel], rot=rot3,
                      y = unit(0, "npc") - tick.units - label.units -
                      labels2.units - unit(6, "mm"),
                      just = c("center", "center"),
                      gp=do.call("gpar", trellis.par.get()$par.xlab.text))
        }
    }

  if (side == "right")
    {
      if (column.panel == cols.per.page)
        {
          axis.options <- axis.options$right
          at2     <- axis.options$at2
          labels2 <- axis.options$labels2
          rot2    <- axis.options$rot2
          ticks2  <- axis.options$ticks2
          at3     <- axis.options$at3
          labels3 <- axis.options$labels3
          rot3    <- axis.options$rot3
          ticks3  <- axis.options$ticks3
          if (is.null(rot2)) rot2 <- 0
          if (is.null(rot3)) rot3 <- 0
          if (is.null(ticks2)) ticks2 <- !is.null(labels2)
          if (is.null(ticks3)) ticks3 <- FALSE

          axis.units <- lattice.getOption("axis.units")$outer$right

          ## space occupied by standard set of ticks
          tick.units <-
            ## if (is.null(scales$at)) unit(0, "lines")
            ## else
              with(axis.units$tick, unit(scales$tck * x, units))

          ## space occupied by standard set of labels, assume 1.2
          ## lines, but we have enough information to compute
          ## explicitly (see how 'lab.grob' is computed in layout.R)

          label.units <-
            if (is.null(scales$at) || !scales$at) unit(0, "lines")
            else unit(1.2, "lines")

          ## want second set of ticks to extend beyond this.

          if (!is.null(at2))
            at2i <-
              if (is.list(at2))
                at2[[row.panel]]
              else at2

          if (!is.null(labels2))
            labels2i <-
              if (is.list(labels2))
                labels2[[row.panel]]
              else labels2

          labels2.units <-
            if (## !is.null(at2) &&  ## Not the same as bottom
                !is.null(labels2))
              unit(1, "grobwidth",
                   data = list(textGrob(unlist(labels2)[1],
                     rot = rot2,
                     y=seq(along=unlist(labels2)[1]),
                     just = if (rot2==90)
                     c("center", "center")
                     else c("left", "center"),
                     gp=do.call("gpar", trellis.par.get()$axis.text))))
            else
              unit(0, "lines")
          labels2.units <- convertWidth(labels2.units, "mm")

          ## recover()
          right.margin <- unit(current.panel.limits()$xlim[2], "native")
          if (ticks2 && (length(at2) > 0) && !is.na(at2))  ## not the same
            grid.segments(y0 = unit(at2i, "native"),
                          y1 = unit(at2i, "native"),
                          x0 = right.margin,
                          x1 = right.margin + tick.units +
                          label.units)

          if ((length(at2) > 0) && !is.na(at2) && length(labels2) > 0)  ## not the same
            grid.text(labels2i, rot = rot2,
                      y = unit(at2i, "native"),
                      x = right.margin + tick.units + label.units +
                      unit(1, "mm"),
                      just = if (rot2==90) c("center", "center") else
                      c("left", "center"),
                      gp=do.call("gpar", trellis.par.get()$axis.text))

          ## other set of labels (without ticks)

          if (ticks3 && !is.null(at3))
            grid.segments(y0 = unit(at3, "native"),
                          y1 = unit(at3, "native"),
                          x0 = right.margin,
                          x1 = right.margin + tick.units +
                          label.units + unit(1, "mm"))

          if (!is.null(at3) && !is.null(labels3))
            grid.text(labels3, rot=rot3,
                      y = unit(at3, "native"),
                      x = right.margin + tick.units +
                      label.units +
                      labels2.units +
                      unit(5, "mm"),
                      just = c("left", "top"),
                      gp=do.call("gpar", trellis.par.get()$par.ylab.text))

          if (is.null(at3) && !is.null(labels3))
            grid.text(labels3[row.panel], rot=rot3,
                      x = right.margin + tick.units +
                      label.units +
                      labels2.units +
                      unit(5, "mm"),
                      just = c("left", "top"),
                      gp=do.call("gpar", trellis.par.get()$par.ylab.text))
        }
    }
}


legendGrob2wt <-   function(...) ## ...is key1, key2, etc
{
  keys <- rev(list(...))

  ## rev is needed because grid and lattice count rows differently (but
  ## only when as.table = FALSE -- you will need to deal with this
  ## somehow, possibly by capturing as.table in your wrapper)

  nkeys <- length(keys)
  key.widths <-
    lapply(keys,
           function(key) unit(1, "grobwidth", data = list(key)))
  key.layout <-
    grid.layout(nrow = nkeys, ncol = 1,
                heights = unit(1, "null"),
                widths = do.call(max, key.widths),
                respect = TRUE)
  key.gf <- frameGrob(layout = key.layout)
  for (i in seq(length=nkeys))
    {
      key.gf <- placeGrob(key.gf, keys[[i]], row = i, col = 1)
    }
  key.gf
}

## source("~/HH-R.package/HH/R/interaction2wt.R")
## source("~/HH-R.package/HH/R/panel.interaction2wt.R")
## trace(axis.i2wt, exit=recover)
