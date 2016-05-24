"panel.cartesian" <-
function(x, y,
         x.label=unique(panel.labels[,"x"]),
         y.label=unique(panel.labels[,"y"]),
         group.label.side="",
         axis3.line=1,
         xg.label, yg.label, g.cex=.7,
         rescale=list(x=TRUE,y=TRUE), ...,
         browser.on=FALSE) {

  lin <- function(x, px, pu.x)
    ((x-px[1])/(px[2]-px[1])) * (pu.x[2]-pu.x[1]) + pu.x[1]

  if.R(r={
    which.parent <- 1
    while(!(exists("rows.per.page", frame=which.parent)))
      which.parent <- which.parent + 1

    cell <- panel.number()

    trellis.object <- get("x", pos=sys.frame(which.parent))
    cols.per.page <- get("cols.per.page", pos=sys.frame(which.parent))
    rows.per.page <- get("rows.per.page", pos=sys.frame(which.parent))
    num.cell <- cols.per.page * rows.per.page
    which.cell <- seq(num.cell)
    this.cell <- match(cell, which.cell)
    if (missing(x.label)) x.label <- ""
    if (missing(y.label)) y.label <- ""
    panel.labels <- cbind(x=rep(x.label, length=cols.per.page),
                          y=rep(rep(y.label, length=rows.per.page),
                            rep(cols.per.page, rows.per.page)))
    these.labels <- panel.labels[this.cell,]
## browser()

    x.up <- trellis.object$x
    y.up <- trellis.object$y

    par.x <- trellis.object$x.scales
    par.y <- trellis.object$y.scales

    pu.x <-
      if (par.x$relation=="same")
        trellis.object$x.limits
      else
        trellis.object$x.limits[[cell]]
    pu.y <-
      if (par.y$relation=="same")
        trellis.object$y.limits
    else
        trellis.object$y.limits[[cell]]
  },
       s={
         ## which.parent <- 1
         ## while(!(exists("cell", frame=sys.parent(which.parent))))
         ##   which.parent <- which.parent + 1

         ## cell <- get("cell", frame=sys.parent(which.parent))
         ## num.cell <- get("num.cell", frame=sys.parent(which.parent))
         ## which.cell <- get("which.cell", frame=sys.parent(which.parent))
         ## this.cell <- match(cell, which.cell)
         ## panel.labels <- get("panel.labels", frame=sys.parent(which.parent))
         ## if (ncol(panel.labels)==1) {
         ##   rn <- dimnames(panel.labels)[[1]]
         ##   panel.labels <-
         ##     if (version$major < 8)
         ##       t(sapply(panel.labels, function(x) sapply(formula(x)[3:2], deparse)))
         ##     else
         ##       do.call("rbind", strsplit(panel.labels[,1], " * "))[,c(3,1)]
         ##   dimnames(panel.labels) <- list(rn, c("x","y"))
         ## }
         ## these.labels <- panel.labels[this.cell,]

         ## x.up <- get("x", frame=sys.parent(which.parent))
         ## y.up <- get("y", frame=sys.parent(which.parent))

         ## glist <- get("glist", frame=sys.parent(which.parent))

         ## par.x <- get("scale.x", frame=sys.parent(which.parent))
         ## par.y <- get("scale.y", frame=sys.parent(which.parent))

         ## gx <- glist[[names(these.labels)[1]]]
         ## rx <- range(x.up[gx$levels[gx$given] == these.labels[1]], na.rm=TRUE)
         ## gy <- glist[[names(these.labels)[2]]]
         ## ry <- range(y.up[gy$levels[gy$given] == these.labels[2]], na.rm=TRUE)

         ## if (any(is.na(rx))) rx <- range(x, na.rm=TRUE)
         ## if (any(is.na(ry))) ry <- range(y, na.rm=TRUE)
       })

##cat(cell,num.cell,which.cell,gy,"\n")
##  if (is.null(gy)) {ry <- c(0,1); y.label <- "abcd"; panel.labels <- cbind(panel.labels,panel.labels)} else

  ##rbind(x.up,y.up,glist$a$given,glist$b$given)

  if.R(r={
    rx <- range(x, na.rm=TRUE)
    if (diff(rx)==0) rx <- rx + c(-1,1)
    pretty.x <- pretty(rx)
    x.cell <- x
    if (rescale$x) {
      px <- pretty.x[c(1, length(pretty.x))]
      px <- px + c(-1,1)*.04*diff(px)
      x.cell <- lin(x, px, pu.x)
      pretty.x.cell <- lin(pretty.x, px, pu.x)
    }
    else {
      x.cell <- x
      pretty.x.cell <- pretty(x)
    }

    ry <- range(y, na.rm=TRUE)
    if (diff(ry)==0) ry <- ry + c(-1,1)
    pretty.y <- pretty(ry)
    y.cell <- y
    if (rescale$y) {
      py <- pretty.y[c(1, length(pretty.y))]
      py <- py + c(-1,1)*.04*diff(py)
      y.cell <- lin(y, py, pu.y)
      pretty.y.cell <- lin(pretty.y, py, pu.y)
    }
    else {
      y.cell <- y
      pretty.y.cell <- pretty(y)
    }
  },
       s={
         if (diff(rx)==0) rx <- rx + c(-1,1)
         pretty.x <- pretty(rx)
         if (rescale$x) {
           px <- pretty.x[c(1, length(pretty.x))]
           px <- px + c(-1,1)*.04*diff(px)
           pu.x <- par()$usr[1:2]
           x.cell <- lin(x, px, pu.x)
           pretty.x.cell <- lin(pretty.x, px, pu.x)
         }
         else {
           x.cell <- x
           pretty.x.cell <- pretty.x
         }

         if (diff(ry)==0) ry <- ry + c(-1,1)
         pretty.y <- pretty(ry)
         if (rescale$y) {
           py <- pretty.y[c(1, length(pretty.y))]
           py <- py + c(-1,1) * .1 * diff(py)
           pu.y <- par()$usr[3:4]
           y.cell <- lin(y, py, pu.y)
           pretty.y.cell <- lin(pretty.y, py, pu.y)
         }
         else {
           y.cell <- y
           pretty.y.cell <- pretty.y
         }
       })

  panel.xyplot(x.cell, y.cell, ...)
##browser()


  if (these.labels[1] == panel.labels[1,1])  ## left column
    if.R(r={
         push.vp.hh()
         if (browser.on) browser()
      grid.yaxis.hh(at=pretty.y.cell, main=TRUE, label=TRUE, labels=pretty.y,
           draw=TRUE)
         pop.vp.hh()
       },
         s=axis(2, at=pretty.y.cell, labels=pretty.y, ticks=TRUE, cex=par.y$cex,
           tck=3*par()$tck, adj=1))

  if (these.labels[1] == rev(unique(panel.labels[,"x"]))[1]) { ## right column
    if.R(r={
         push.vp.hh()
         if (browser.on) browser()
      grid.yaxis.hh(at=pretty.y.cell, main=FALSE, label=FALSE, draw=TRUE)
      grid.text(x=unit(pu.x[2],"native"),
                y=unit(mean(pu.y),"native"),
                just="left", draw=TRUE,
                paste("    ",
                as.character(y.label[match(these.labels[2],
                                           unique(panel.labels[,2]))])))
         pop.vp.hh()
    },
         s={
           axis(4, at=pretty.y.cell, labels=FALSE, ticks=TRUE, cex=par.y$cex,
                tck=3*par()$tck)
           mtext(as.character(y.label[match(these.labels[2],
                                           unique(panel.labels[,2]))]),
                side=4, line=3, at=mean(par()$usr[3:4]), srt=0, adj=1, cex=g.cex)
  })}

  if (these.labels[2] == panel.labels[1,2]) ## bottom row
    if.R(r={
         push.vp.hh()
         if (browser.on) browser()
         grid.xaxis.hh(at=pretty.x.cell, main=TRUE, label=TRUE, labels=pretty.x, draw=TRUE)
         pop.vp.hh()
       },
         s=
    axis(1, at=pretty.x.cell, labels=pretty.x, ticks=TRUE, cex=par.x$cex,
         tck=3*par()$tck))

  if (these.labels[2] == rev(unique(panel.labels[,"y"]))[1]) { ## top row
    if.R(r={
         push.vp.hh()
         if (browser.on) browser()
           c.v <- current.viewport()
           c.v$height <- unit(1,"npc")+unit(2*axis3.line,"lines")
           c.v$width <- unit(1,"npc")
           pushViewport(c.v)
           grid.xaxis.hh(at=pretty.x.cell, main=FALSE, label=FALSE)
         if (browser.on) browser()
           popViewport()
         grid.text(as.character(x.label[match(these.labels[1],
                                           unique(panel.labels[,1]))]),
                y=unit(pu.y[2],"native")+unit(3.7*axis3.line,"lines"),
                x=unit(mean(pu.x),"native"), draw=TRUE)
      if ((match("group", names(these.labels), 0)) && ## group exists
          (these.labels[1] == panel.labels[1,1]))     ## upper left corner
        if (group.label.side=="left")
          grid.text(these.labels["group"],
                    x=unit(3, "lines"), y=unit(pu.y[1],"native"),
                    gp=gpar(outer=TRUE, adj=1, srt=90))
        else ## (group.label.side=="top")
          grid.text(these.labels["group"],
                    x=unit(pu.x[2],"native"), y=unit(1.5, "lines"),
                    gp=gpar(adj=0, srt=0))
         pop.vp.hh()
       },
         s={
           axis(3, at=pretty.x.cell, labels=FALSE, ticks=TRUE, cex=par.x$cex,
                line=axis3.line, tck=3*par()$tck)
           mtext(as.character(x.label[match(these.labels[1],
                                            unique(panel.labels[,1]))]),
                 side=3, line=axis3.line+.5, at=mean(par()$usr[1:2]), cex=g.cex)
           if ((match("group", names(these.labels), 0)) && ## group exists
               (these.labels[1] == panel.labels[1,1]))     ## upper left corner
             if (group.label.side=="") {}
             else
               if (group.label.side=="left")
                 mtext(these.labels["group"],
                       side=2, outer=TRUE, line=3, adj=1, srt=90, at=par()$usr[3])
               else ## (group.label.side=="top")
                 mtext(these.labels["group"],
                       side=3, outer=FALSE, line=1.5, adj=0, srt=0, at=par()$usr[2])
         })
  }

  if (cell==1) {
    if.R(
         r={
         push.vp.hh()
         if (browser.on) browser()
           if (!missing(xg.label))
             grid.text(xg.label,
                       x=unit(pu.x[1],"native"), y=unit(pu.y[2],"native"),
                       gp=gpar(outer=TRUE))
           if (!missing(yg.label))
            grid.text(yg.label,
                       x=unit(pu.x[2],"native"), y=unit(pu.y[1],"native"),
                       gp=gpar(outer=TRUE))
         pop.vp.hh()
         },
         s={
           if (!missing(xg.label))
             mtext(xg.label, side=3, line=axis3.line, outer=TRUE)
           if (!missing(yg.label))
             mtext(yg.label, side=4, outer=TRUE)
         }
         )
  }
}
