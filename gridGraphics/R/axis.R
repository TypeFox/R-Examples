
# C_axis(side, at, labels, 
#        tick, line, pos, outer, font, lty, lwd, lwd.ticks, col, 
#        col.ticks, hadj, padj, ...)
C_axis <- function(x) {
    dev.set(recordDev())
    # Blank out x$cex because we want par$cex to be par$cexbase
    x$cex <- NULL
    par <- currentPar(x[-(1:16)])
    dev.set(playDev())
    depth <- gotovp(NA)
    side <- x[[2]]
    if (is.null(x[[3]])) {
        ticks <- defaultTicks(side, par)
    } else {
        ticks <- x[[3]]
    }
    labels <- x[[4]]
    doticks <- x[[5]]
    if (is.na(doticks))
        doticks <- TRUE
    # [1] to guard against 'line' vector of length > 1
    line <- x[[6]][1]
    pos <- x[[7]]
    # Not sure on the logic of the following, just emulating C code
    lineoff <- 0
    if (!is.finite(line)) {
        line <- par$mgp[3]
        lineoff <- line
    } 
    if (is.finite(pos))
        lineoff <- 0
    outer <- x[[8]]
    font <- FixupFont(x[[9]], NA)
    lty <- FixupLty(x[[10]], 0)
    lwd <- FixupLwd(x[[11]], 1)
    lwd.ticks <- FixupLwd(x[[12]], 1)
    col <- FixupCol(x[[13]], par$fg, par$bg)
    col.ticks <- FixupCol(x[[14]], col, par$bg)
    hadj <- x[[15]]
    padj <- x[[16]]
    # NOTE: the use of 'trim=TRUE' in format() to mimic use of,
    #       e.g., EncodeReal0(), within labelformat() in plot.c
    if (is.null(labels)) {
        drawLabels <- TRUE
        labels <- format(ticks, trim=TRUE)
    } else if (is.logical(labels)) {
        if (labels) {
            drawLabels <- TRUE
            labels <- format(ticks, trim=TRUE)
        } else {
            drawLabels <- FALSE
        }
    } else {
        drawLabels <- TRUE
    }
    if (is.finite(par$tck)) {
        if (par$tck > 0.5) {
            tickLength <- unit(par$tck, "npc")
        } else {
            tickLength <- min(convertWidth(unit(par$tck, "npc"), "in"),
                              convertHeight(unit(par$tck, "npc"), "in"))
        }
    } else {
        tickLength <- unit(par$cin[2]*par$tcl*par$cex, "in")
    }
    returnvp <- NULL
    if (side == 1 && par$xaxt != "n") {
        if (is.finite(pos)) {
            axis_base <- unit(pos, "native")
        } else {
            if (outer) {
                returnvp <- pushTempViewport(side)
            } 
            axis_base <- unit(0, "npc") - unit(line*par$cin[2]*par$cex, "in")
        }
        # Now that have generated tick labels from tick locations
        # (if necessary), can transform tick locations (if necessary)
        # for log transforms
        ticks <- tx(ticks, par)
        # Clip ticks (and labels) to plot boundaries
        ticksub <- ticks >= par$usr[1] & ticks <= par$usr[2]
        if (doticks) {
            if (lwd > 0) {
                grid.segments(unit(min(par$usr[2], max(par$usr[1], min(ticks))),
                                   "native"),
                              axis_base,
                              unit(min(par$usr[2], max(par$usr[1], max(ticks))),
                                   "native"),
                              axis_base,
                              gp=gpar(col=col, lwd=lwd, lty=lty,
                                  lineend=par$lend, linemitre=par$lmitre,
                                  linejoin=par$ljoin),
                              name=grobname("bottom-axis-line"))
            }
            if (lwd.ticks > 0) {
                grid.segments(unit(ticks[ticksub], "native"),
                              axis_base,
                              unit(ticks[ticksub], "native"),
                              axis_base + tickLength,
                              gp=gpar(col=col.ticks, lwd=lwd.ticks, lty=lty,
                                  lineend=par$lend, linemitre=par$lmitre,
                                  linejoin=par$ljoin),
                              name=grobname("bottom-axis-ticks"))
            }
        }
        # NOTE: the following includes calculation based on par(mgp)
        #       to get the margin line to draw on
        #       PLUS adjustment made in GMtext() based on that line value
        if (drawLabels) {
            labLine <- - (convertY(axis_base, "in", valueOnly=TRUE)/
                          (par$cin[2]*par$cex)) +
                       par$mgp[2] - lineoff
            GMtext(labels[ticksub], 1, line=labLine,
                   at=unit(ticks[ticksub], "native"), las=par$las, 
                   xadj=computeXAdj(hadj, side, par$las),
                   yadj=computePAdj(padj, side, par$las),
                   mex=par$mex, cin=par$cin, 
                   cex=par$cex.axis*par$cex, linecex=par$mex*par$cex,
                   font=par$font.axis, family=par$family,
                   col=par$col.axis, lheight=par$lheight,
                   yLineBias=par$ylbias, allowOverlap=FALSE,
                   label="bottom-axis-labels")
        }
    } else if (side == 2 && par$yaxt != "n") {
        if (is.finite(pos)) {
            axis_base <- unit(pos, "native")
        } else {
            if (outer) {
                returnvp <- pushTempViewport(side)
            } 
            axis_base <- unit(0, "npc") - unit(line*par$cin[2]*par$cex, "in")
        }
        ticks <- ty(ticks, par)
        ticksub <- ticks >= par$usr[3] & ticks <= par$usr[4]
        if (doticks) {
            if (lwd > 0) {
                grid.segments(axis_base,
                              unit(min(par$usr[4], max(par$usr[3], min(ticks))),
                                   "native"),
                              axis_base,
                              unit(min(par$usr[4], max(par$usr[3], max(ticks))),
                                   "native"),
                              gp=gpar(col=col, lwd=lwd, lty=lty,
                                  lineend=par$lend, linemitre=par$lmitre,
                                  linejoin=par$ljoin),
                              name=grobname("left-axis-line"))
            }
            if (lwd.ticks > 0) {            
                grid.segments(axis_base,
                              unit(ticks[ticksub], "native"),
                              axis_base + tickLength,
                              unit(ticks[ticksub], "native"),
                              gp=gpar(col=col.ticks, lwd=lwd.ticks, lty=lty,
                                  lineend=par$lend, linemitre=par$lmitre,
                                  linejoin=par$ljoin),
                              name=grobname("left-axis-ticks"))
            }
        }
        if (drawLabels) {
            labLine <- - (convertX(axis_base, "in", valueOnly=TRUE)/
                          (par$cin[2]*par$cex)) +
                       par$mgp[2] - lineoff
            GMtext(labels[ticksub], 2, line=labLine,
                   at=unit(ticks[ticksub], "native"), las=par$las,
                   xadj=computeXAdj(hadj, side, par$las),
                   yadj=computePAdj(padj, side, par$las),
                   mex=par$mex, cin=par$cin,
                   cex=par$cex.axis*par$cex, linecex=par$mex*par$cex,
                   font=par$font.axis, family=par$family,
                   col=par$col.axis, lheight=par$lheight,
                   yLineBias=par$ylbias, allowOverlap=FALSE,
                   label="left-axis-labels")
        }
    } else if (side == 3 && par$xaxt != "n") {
        if (is.finite(pos)) {
            axis_base <- unit(pos, "native")
        } else {
            if (outer) {
                returnvp <- pushTempViewport(side)
            } 
            axis_base <- unit(1, "npc") + unit(line*par$cin[2]*par$cex, "in")
        }
        ticks <- tx(ticks, par)
        ticksub <- ticks >= par$usr[1] & ticks <= par$usr[2]
        if (doticks) {
            if (lwd > 0) {
                grid.segments(unit(min(par$usr[2], max(par$usr[1], min(ticks))),
                                   "native"),
                              axis_base,
                              unit(min(par$usr[2], max(par$usr[1], max(ticks))),
                                   "native"),
                              axis_base,
                              gp=gpar(col=col, lwd=lwd, lty=lty,
                                  lineend=par$lend, linemitre=par$lmitre,
                                  linejoin=par$ljoin),
                              name=grobname("top-axis-line"))
            }
            if (lwd.ticks > 0) {
                grid.segments(unit(ticks[ticksub], "native"),
                              axis_base,
                              unit(ticks[ticksub], "native"),
                              axis_base - tickLength,
                              gp=gpar(col=col.ticks, lwd=lwd.ticks, lty=lty,
                                  lineend=par$lend, linemitre=par$lmitre,
                                  linejoin=par$ljoin),
                              name=grobname("top-axis-ticks"))
            }
        }
        if (drawLabels) {
            labLine <- (convertY(axis_base, "in", valueOnly=TRUE)/
                        (par$cin[2]*par$cex)) +
                       par$mgp[2] - lineoff -
                       (convertY(unit(1, "npc"), "in", valueOnly=TRUE)/
                        (par$cin[2]*par$cex))
            GMtext(labels[ticksub], 3, line=labLine,
                   at=unit(ticks[ticksub], "native"), las=par$las, 
                   xadj=computeXAdj(hadj, side, par$las),
                   yadj=computePAdj(padj, side, par$las),
                   mex=par$mex, cin=par$cin, 
                   cex=par$cex.axis*par$cex, linecex=par$mex*par$cex,
                   font=par$font.axis, family=par$family,
                   col=par$col.axis, lheight=par$lheight,
                   yLineBias=par$ylbias, allowOverlap=FALSE,
                   label="top-axis-labels")
        }
    } else if (side == 4 && par$yaxt != "n") {
        if (is.finite(pos)) {
            axis_base <- unit(pos, "native")
        } else {
            if (outer) {
                returnvp <- pushTempViewport(side)
            } 
            axis_base <- unit(1, "npc") + unit(line*par$cin[2]*par$cex, "in")
        }
        ticks <- ty(ticks, par)
        ticksub <- ticks >= par$usr[3] & ticks <= par$usr[4]
        if (doticks) {
            if (lwd > 0) {
                grid.segments(axis_base,
                              unit(min(par$usr[4], max(par$usr[3], min(ticks))),
                                   "native"),
                              axis_base,
                              unit(min(par$usr[4], max(par$usr[3], max(ticks))),
                                   "native"),
                              gp=gpar(col=col, lwd=lwd, lty=lty,
                                  lineend=par$lend, linemitre=par$lmitre,
                                  linejoin=par$ljoin),
                              name=grobname("right-axis-line"))
            }
            if (lwd.ticks > 0) {
                grid.segments(axis_base,
                              unit(ticks[ticksub], "native"),
                              axis_base - tickLength,
                              unit(ticks[ticksub], "native"),
                              gp=gpar(col=col.ticks, lwd=lwd.ticks, lty=lty,
                                  lineend=par$lend, linemitre=par$lmitre,
                                  linejoin=par$ljoin),
                              name=grobname("right-axis-ticks"))
            }
        }
        if (drawLabels) {
            labLine <- (convertX(axis_base, "in", valueOnly=TRUE)/
                        (par$cin[2]*par$cex)) +
                       par$mgp[2] - lineoff -
                       (convertX(unit(1, "npc"), "in", valueOnly=TRUE)/
                        (par$cin[2]*par$cex))
            GMtext(labels[ticksub], 4, line=labLine,
                   at=unit(ticks[ticksub], "native"), las=par$las,
                   xadj=computeXAdj(hadj, side, par$las),
                   yadj=computePAdj(padj, side, par$las),
                   mex=par$mex, cin=par$cin,
                   cex=par$cex.axis*par$cex, linecex=par$mex*par$cex,
                   font=par$font.axis, family=par$family,
                   col=par$col.axis, lheight=par$lheight,
                   yLineBias=par$ylbias, allowOverlap=FALSE,
                   label="right-axis-labels")
        }
    }
    # Undo any temporary "outer" viewport
    if (!is.null(returnvp)) {
        upViewport()
        downViewport(returnvp)
    }
    upViewport(depth)
}

defaultTicks <- function(side, par) {
    axp <- switch(side, par$xaxp, par$yaxp, par$xaxp, par$yaxp)
    usr <- switch(side, par$usr[1:2], par$usr[3:4], par$usr[1:2], par$usr[3:4])
    log <- switch(side, par$xlog, par$ylog, par$xlog, par$ylog)
    axTicks(side, axp, usr, log)
}

computeXAdj <- function(hadj, side, las) {
    if (is.finite(hadj)) {
        xadj <- hadj
    } else {
        if (side == 1 || side == 3) {
            if (las == 2 || las == 3) {
                if (side == 1) {
                    xadj <- 1
                } else {
                    xadj <- 0
                }
            } else {
                xadj <- 0.5
            }
        } else {
            if (las == 1 || las == 2) {
                if (side == 2) {
                    xadj <- 1
                } else {
                    xadj <- 0
                }
            } else {
                xadj <- 0.5
            }
        } 
    }
    xadj
}

computePAdj <- function(padj, side, las) {
    if (!is.finite(padj)) {
        padj <- switch(las + 1,
                       0,
                       switch(side, 0, 0.5, 0, 0.5),
                       0.5,
                       switch(side, 0.5, 0, 0.5, 0))
    }
    padj
}

pushTempViewport <- function(side) {
    cvp <- current.viewport()
    xscale <- cvp$xscale
    yscale <- cvp$yscale
    name <- cvp$name
    returnvp <- upViewport(2)
    if (side == 1 || side == 3) {
        pushViewport(viewport(x=grobX(rectGrob(vp=returnvp), "west"),
                              width=grobWidth(rectGrob(vp=returnvp)),
                              just="left",
                              xscale=xscale,
                              name=paste0(name, "-outer-axis")))
    } else { # side == 2 || side == 4
        pushViewport(viewport(y=grobY(rectGrob(vp=returnvp), "south"),
                              height=grobHeight(rectGrob(vp=returnvp)),
                              just="bottom",
                              yscale=yscale,
                              name=paste0(name, "-outer-axis")))        
    }
    returnvp
}
