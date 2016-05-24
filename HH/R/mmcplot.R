mmcplot <- function (mmc, ...)
  UseMethod("mmcplot")


mmcplot.mmc <- function(mmc, col=col, lwd=lwd, lty=lty, ...,
                        style=c("isomeans", "confint", "both"),
                        type=c("mca", "lmat", "linfct", "none")) {
  style <- match.arg(style)
  type <- match.arg(type)
  if (type=="linfct")
    type <- "lmat"
  if (type %in% c("mca", "lmat"))
    switch(style,
           isomeans=mmcisomeans(mmc, col=col, lwd=lwd, lty=lty, type=type, ...),
           confint=mmcmatch(mmc, col=col, lwd=lwd, lty=lty, type=type, ...),
           both=mmcboth(mmc, col=col, lwd=lwd, lty=lty, type=type, ...))
  else
    ## mmcglht(mmc[[type]], ...)
    mmcmatch(mmc, type=type, col=col, lwd=lwd, lty=lty, ...)
}

mmcplot.glht <- function(mmc, col=c("black","red"), lwd=c(1,1), lty=c(2,1), focus=mmc$focus, ...) {
  if (is.null(focus))
    stop("Please specify 'focus='.", call.=FALSE)
  mmcmatch(as.multicomp.glht(mmc, focus=focus, ...), col=col, lwd=lwd, lty=lty, ..., xlim.match=FALSE)
}

mmcplot.mmc.multicomp <- function(mmc, col=c("black","red"), lwd=c(1,1), lty=c(2,1), ...)
  mmcplot.mmc(mmc, col=col, lwd=lwd, lty=lty, ..., mmc.object.name=deparse(substitute(mmc)))

mmcplot.multicomp <- function(mmc, col=col, lwd=lwd, lty=lty, ...)
  mmcmatch(mmc, col=col, lwd=lwd, lty=lty, ..., xlim.match=FALSE)

## mmcplot.multicomp.hh <- function(mmc, ...)
##   mmcmatch(mmc, ..., xlim.match=FALSE)

mmcplot.default <- function(mmc, ...)
  stop("Only mmc and glht objects can be plotted by mmcplot.", call. = FALSE)



mmcisomeans <- function(mmc, col=col, lwd=lwd, lty=lty, type="mca", xlim=NULL, ylim=NULL, ...,
                        axis.right=2.2,
                        ylab=paste(
                          mmc$none$ylabel, "means",
                          " | ",
                          mmc$none$focus, "level"),
                        ylab.right=NULL,
                        xlab="contrast value",
                        contrast.label=TRUE,
                        means.height=TRUE) {
  mmc.type <- mmc[[type]]
  if (is.null(mmc.type)) {
    if (type %in% c("lmat", "linfct"))
      stop("type='", type, "' not in ", list(...)$mmc.object.name, call.=FALSE)
    else
      stop("invalid type", call.=FALSE)
  }
  tmp <- data.frame(mmc.type$table, height=mmc.type$height, height.2=mmc.type$height/2,
                    contrast.name=rownames(mmc.type$table), row=1:length(mmc.type$height),
                    signif=factor(
                      (mmc.type$table[,"lower"] * mmc.type$table[,"upper"] > 0),
                      levels=c(FALSE, TRUE))
                    )
  tmp

  max.upper <- max(sapply(mmc[-2],
                          function(mmci) max(mmci$table[,"upper"]))) ## "none" must be in position 2
  min.lower <- min(sapply(mmc[-2],
                          function(mmci) min(mmci$table[,"lower"]))) ## "none" must be in position 2

  means <- mmc$none$table[,"estimate"]

  if (is.null(ylim)) {
    ylim <- range(means)
    ylim <- ylim + diff(ylim)*c(-.15, .05)
  }

  if (is.infinite(max.upper)) max.upper <-  diff(ylim)
  if (is.infinite(min.lower)) min.lower <- -diff(ylim)
  if (is.null(xlim))
    xlim <- range(c(min.lower, -diff(ylim), diff(ylim), max.upper))*1.05

  aspect.ratio <- 2*diff(ylim)/diff(xlim)


  xyplot(height.2 ~ estimate, groups=signif, data=tmp,
         means=means,
         lower=tmp$lower,
         upper=tmp$upper,
         contrast.name=tmp$contrast.name,
         par.settings=list(clip=list(panel=FALSE), layout.widths=list(axis.right=axis.right)),
         col=col, lwd=lwd, lty=lty,
         scales=list(y=list(
                       at=if (means.height) means else NULL,
                       tck=c(TRUE, FALSE))),
         xlim=xlim,
         ylim=ylim,
         aspect=aspect.ratio,
         ylab=ylab,
         ylab.right=ylab.right,
         xlab=xlab,
         panel=function(..., ylabel, focus) {
           panel.isomeans(ybar=means, ...)
           panel.superpose(...)
           panel.axis("left", at=means, labels=names(means), outside=FALSE, half=FALSE)
           if (contrast.label)
             panel.axis("bottom", at=current.panel.limits()$xlim[2],
                        labels="                       contrasts",
                        ## extra space because adjust not available
                        rot=0, outside=TRUE, tck=FALSE)
         },
         panel.groups=panel.confintMMC,
         ...)
}


mmcmatch <- function(mmc, col=col, lwd=lwd, lty=lty, type="mca", xlim=NULL, ylim=NULL, ...,
                     axis.right=2.2,
                     ylab=NULL,
                     ylab.right=NULL,
                     xlab="contrast value",
                     contrast.label=TRUE,
                     xlim.match=(type != "none")) {
  mmc.type <- mmc[[type]]
  if (is.null(mmc.type)) mmc.type <- mmc
  if (is.null(mmc.type)) {
    if (type %in% c("lmat", "linfct"))
      stop("type='", type, "' not in ", list(...)$mmc.object.name, call.=FALSE)
    else
      stop("invalid type", call.=FALSE)
  }

  tmp <- data.frame(mmc.type$table, height=mmc.type$height, height.2=mmc.type$height/2,
                    contrast.name=rownames(mmc.type$table), row=1:length(mmc.type$height),
                    signif=factor(
                      (mmc.type$table[,"lower"] * mmc.type$table[,"upper"] > 0),
                      levels=c(FALSE, TRUE))
                    )
  tmp

  if (names(mmc)[[2]] == "none") {
    max.upper <- max(sapply(mmc[-2],
                            function(mmci) max(mmci$table[,"upper"]))) ## "none" must be in position 2
    min.lower <- min(sapply(mmc[-2],
                            function(mmci) min(mmci$table[,"lower"]))) ## "none" must be in position 2
  }
  else {
    max.upper <- nrow(tmp)
    min.lower <- 1
  }

  means <- mmc$none$table[,"estimate"]
  if (is.null(means)) means <- tmp$estimate
  if (is.null(ylim)) {
    ylim <- range(means)
    ylim <- ylim + diff(ylim)*c(-.05, .05)
  }

  if (is.infinite(max.upper)) max.upper <-  diff(ylim)
  if (is.infinite(min.lower)) min.lower <- -diff(ylim)
  if (is.null(xlim) && xlim.match)
    xlim <- range(c(min.lower, -diff(ylim), diff(ylim), max.upper))*1.05
  if (is.null(xlim) && !xlim.match) {
    xlim <- range(tmp$lower, tmp$upper)
    if (is.infinite(xlim[2])) xlim[2] <- max(tmp$estimate + (tmp$estimate-tmp$lower))
    if (is.infinite(xlim[1])) xlim[1] <- max(tmp$estimate - (tmp$upper-tmp$estimate))
    xlim <- xlim + diff(xlim)*c(-.05, .05)
  }

  aspect.ratio <- 2*diff(ylim)/diff(xlim)


  xyplot(rev(row) ~ estimate, groups=signif, data=tmp,
         means=means,
         lower=tmp$lower,
         upper=tmp$upper,
         contrast.name=tmp$contrast.name,
         par.settings=list(clip=list(panel=FALSE), layout.widths=list(axis.right=axis.right)),
         col=col, lwd=lwd, lty=lty,
         ## scales=list(y=list(alternating=0, at=NULL)),
         scales=list(y=list(alternating=ifelse(xlim.match, 1, 0),
                       at=1:length(tmp$height.2), label=rev(signif(tmp$height.2, 4)))),
         xlim=xlim,
         ylim=range(tmp$row) + c(-.5, .5),
         ylab=ylab,
         ylab.right=ylab.right,
         xlab=xlab,
         panel=function(...,
           lty.contr0=2,
           col.contr0='darkgray',
           lwd.contr0=1
           ) {
           panel.abline(v=0, lty=lty.contr0, col=col.contr0, lwd=lwd.contr0)
           if (contrast.label)
             panel.axis("bottom", at=current.panel.limits()$xlim[2],
                        labels="                       contrasts",
                        ## extra space because adjust not available
                        rot=0, outside=TRUE, tck=FALSE)
           if (contrast.label && xlim.match)
             panel.axis("bottom", at=current.panel.limits()$xlim[1],
                        labels="height                       ",
                        ## extra space because adjust not available
                        rot=0, outside=TRUE, tck=FALSE)
           panel.superpose(...)
         },
         panel.groups=panel.confintMMC,
         ...)
}


mmcboth <- function(mmc, col=col, lwd=lwd, lty=lty, type="mca", h=c(.7, .3), xlim=NULL, ylim=NULL, ...,
                    ylab.right=NULL, MMCname="MMC", Tiebreakername="Tiebreaker") {
  aa <- mmcisomeans(mmc=mmc, col=col, lwd=lwd, lty=lty, type=type, xlim=xlim, ylim=ylim, ..., contrast.label=TRUE)
  bb <- mmcmatch(mmc=mmc, col=col, lwd=lwd, lty=lty, type=type, xlim=xlim, ylim=ylim, ..., contrast.label=TRUE, xlim.match=TRUE)
  cc0 <- c(MMC=aa, Tiebreaker=bb, layout=c(1,2))
  if (!missing(MMCname)) cc0$condlevels[[1]][1] <- MMCname
  if (!missing(Tiebreakername)) cc0$condlevels[[1]][2] <- Tiebreakername
  cc <- resizePanels(cc0, h=h)
  dd <- (update(cc, as.table=TRUE, between=list(y=1)))
  dd$aspect.ratio <- mmcAspect(dd) / h[1]

  dd$y.scales$at=list(aa$y.scales$at, bb$y.scales$at)
  dd$y.scales$labels=list(aa$y.scales$labels, bb$y.scales$labels)
  dd$y.scales$rot=0
  dd$ylab <- c(" ", aa$ylab)
  dd$ylab.right <- ylab.right
  if (dd$x.scales$relation != "same") {
    dd$x.limits <- dd$x.limits[[1]]
    dd$x.scales$relation <- "same"
  }
  dd
}


panel.confintMMC <- function(x, y, subscripts, ..., col, lwd, lty, lower, upper,
                             contrast.name, right.text.cex=.8,
                             contrast.height=FALSE) {
  if (is.infinite(lower[subscripts[1]]))
    left <- current.panel.limits()$xlim[1]
  else
    left <- lower[subscripts]
  if (is.infinite(upper[subscripts[1]]))
    right <- current.panel.limits()$xlim[2]
  else
    right <- upper[subscripts]
  panel.segments(left, y, right, y, col=col, lwd=lwd, lty=lty)
  panel.points(x, y, pch=3, col=col)
  panel.axis("right", at=y, labels=contrast.name[subscripts],
             text.col=col, line.col=col, outside=TRUE, half=FALSE, text.cex=right.text.cex)
  if (contrast.height)
      panel.axis("left", at=y,
             text.col=col, line.col=col, outside=TRUE, half=FALSE, text.cex=right.text.cex)

}


panel.isomeans <- function(ybar,
                           lty.iso=7,
                           col.iso='darkgray',
                           lwd.iso=1,
                           lty.contr0=2,
                           col.contr0='darkgray',
                           lwd.contr0=1,
                           ...,
                           col, lwd, lty ## capture potentially ambiguous name
                           ) {
  ## the additional feature of the Hsu Peruggia plot is the
  ## iso-lines for each level of the factor
  ## names(ybar) <- seq(along=ybar)
  ## ybaro <- order(ybar)
  ## ybar <- ybar[ybaro]
  ybar.min <- min(ybar)
  ybar.max <- max(ybar)
  panel.segments(ybar.min-ybar, (ybar.min+ybar)/2,
                 ybar.max-ybar, (ybar.max+ybar)/2, lty=lty.iso, col=col.iso, lwd=lwd.iso)
  panel.segments(ybar-ybar.min, (ybar+ybar.min)/2,
                 ybar-ybar.max, (ybar+ybar.max)/2, lty=lty.iso, col=col.iso, lwd=lwd.iso)
  panel.abline(v=0, lty=lty.contr0, col=col.contr0, lwd=lwd.contr0)

  panel.axis("bottom", at=c(ybar.min-ybar, ybar-ybar.min), labels=FALSE, half=FALSE)

  ydif <- .025*diff(current.panel.limits()$ylim)

  ybar.names <- names(ybar)
  ybar.min.index <- match(ybar.min, ybar)

  panel.text( ybar.min-ybar[-ybar.min.index],
             (ybar+ybar.min)[-ybar.min.index]/2-ydif,
             ybar.names[-ybar.min.index], adj=.9, col=col.iso)
  panel.text(-ybar.min+ybar[-ybar.min.index],
             (ybar+ybar.min)[-ybar.min.index]/2-ydif,
             ybar.names[-ybar.min.index], adj=.1, col=col.iso)
  panel.text(0                              ,
             (ybar+ybar.min)[ ybar.min.index]/2-ydif,
             ybar.names[ ybar.min.index], adj=.5, col=col.iso)
}


mmcAspect <- function(trellis) {
  diff.y <- if (is.list(trellis$y.limits))
    diff(trellis$y.limits[[1]])
  else
    diff(trellis$y.limits)
  diff.x <- if(is.list(trellis$x.limits))
    diff(trellis$x.limits[[1]])
  else
    diff(trellis$x.limits)
  (2 * diff.y / diff.x)
}

## ## still needed
## 1. color arguments:
##    significant, non-significant pairwise, non-significant lmat, isomeans, strip background
## 2. left axis: level names or numbers with legend
## 3. main, sub, xlab, xlab.top, ylab, ylab.right
## 4. right-axis outer=FALSE for lmat contrasts


## this doesn't work
## xyplot(1:10 ~ 1:10,
##          ylab=expression(
##            atop(as.expression(bquote(.(mmc$none$ylabel))), "means") ~
##            "|" ~
##            atop(as.expression(bquote(.(mmc$none$focus)), "level")))
##        )

## idiom works
## xyplot(0 ~ 1, sub=as.expression(
##         bquote(pi==.(pi) ~ e==.(exp(1)))
## ))


mmcPruneIsomeans <- function(mmc, keep=NULL) {
  if (is.null(keep)) return(mmc)

  ## for non-null keep
  mmc$none$table <- mmc$none$table[keep,]
  mmc$none$height <- mmc$none$height[keep]
  mmc
}

