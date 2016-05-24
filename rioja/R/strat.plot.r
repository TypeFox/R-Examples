strat.plot <- function (d, yvar = NULL, scale.percent = FALSE, graph.widths=1, minmax=NULL, scale.minmax = TRUE,
    xLeft = 0.07, xRight = 1, yBottom = 0.07, yTop = 0.8, title = "", cex.title=1.8, y.axis=TRUE,
    min.width = 5, ylim = NULL, y.rev = FALSE, y.tks=NULL, ylabel = "", cex.ylabel=1, cex.yaxis=1,
    xSpace = 0.01, x.pc.inc=10, x.pc.lab=TRUE, x.pc.omit0=TRUE, wa.order = "none", plot.line = TRUE, col.line = "black", lwd.line = 1, 
    plot.bar = TRUE, lwd.bar = 1, col.bar = "grey", sep.bar = FALSE, 
    plot.poly = FALSE, col.poly = "grey", col.poly.line = NA, lwd.poly = 1,
    plot.symb = FALSE, symb.pch=19, symb.cex=1, x.names=NULL,
    cex.xlabel = 1.1, srt.xlabel=90, mgp=NULL, cex.axis=.8, clust = NULL, clust.width=0.1,
    orig.fig=NULL, add.smooth=FALSE, col.smooth="red", smooth.span=0.1, lwd.smooth=1, 
    exag=FALSE, exag.mult=5, col.exag="grey90", exag.alpha=0.2, add=FALSE, ...)
{
    fcall <- match.call(expand.dots=TRUE)
    if (!is.null(clust)) {
     if (class(clust)[1]!="chclust") 
        stop("clust must be a chclust object")
    }
    if (!is.null(clust)) 
        xRight = xRight - clust.width
    if (is.null(yvar)) {
        yvar <- 1:nrow(d)
        if (is.null(ylim)) {
           ylim=c(0.5, nrow(d)+0.5)
        }
    }
    if (is.null(x.names))
       x.names=colnames(d)   
    if (is.null(ylim)) 
        ylim = c(min(yvar, na.rm=TRUE), max(yvar, na.rm=TRUE))
    oldfig = par("fig")
    oldmai <- par("mai")
    if (is.null(orig.fig)) {
       orig.fig = par("fig")
    }
    nsp <- ncol(d)
    nsam <- nrow(d)
    if (scale.percent==TRUE & length(x.pc.inc) > 1) {
       if (length(x.pc.inc) != nsp)
          stop("length of x.pc.inc should equal number of curves")
    } else {
       x.pc.inc <- rep(x.pc.inc[1], nsp)
    }
    if (!is.null(minmax)) {
       if (ncol(minmax) != 2)
          stop("minmax should have 2 columns")
       if (nrow(minmax) != nsp)
          stop("number of rows of minmax should equal number of curves")
    }
    par(mai = c(0, 0, 0, 0))
    if (length(graph.widths) == 1)
       graph.widths <- rep(1, nsp)
    if (length(graph.widths) != nsp)
       stop("Length of graph.widths should equal number of curves")
    if (length(exag) == 1)
      exag <- rep(exag[1], nsp)
    if (length(exag) != nsp)
      stop("Length of exag should equal number of curves")
    if (length(exag.mult) == 1)
      exag.mult <- rep(exag.mult[1], nsp)
    if (length(exag.mult) != nsp)
      stop("Length of exag.mult should equal number of curves")
    if (length(col.exag) == 1)
      col.exag <- rep(col.exag[1], nsp)
    if (length(col.exag) != nsp)
      stop("Length of col.exag should equal number of curves")
    if (length(add.smooth) == 1)
      add.smooth <- rep(add.smooth[1], nsp)
    if (length(add.smooth) != nsp)
      stop("Length of add.smooth should equal number of curves")
    if (length(smooth.span) == 1)
      smooth.span <- rep(smooth.span[1], nsp)
    if (length(smooth.span) != nsp)
      stop("Length of smooth.span should equal number of curves")
    cc.line <- rep(col.line, length.out=nsp)
    if (sep.bar)
      cc.bar <- rep(col.bar, length.out=nsam)
    else      
      cc.bar <- rep(col.bar, length.out=nsp)
    cc.poly <- rep(col.poly, length.out=nsp)
    cc.poly.line <- rep(col.poly.line, length.out=nsp)
    if(plot.poly)
      plot.line <- FALSE
    make.col <- function(x, alpha) {
      apply(col2rgb(x)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha))
    }
    if (col.exag[1] == "auto")
      col.exag <- make.col(cc.poly, exag.alpha)
    inc <- 0.002
    if (wa.order == "topleft" || wa.order == "bottomleft") {
        colsum <- colSums(d)
        opt <- (t(d) %*% yvar)/colsum
        if ((wa.order == "topleft" & !y.rev) | (wa.order == "bottomleft" & y.rev))
            opt.order <- rev(order(opt))
        else opt.order <- order(opt)
        d <- d[, opt.order]
        if (!is.null(minmax)) 
           minmax <- minmax[opt.order, ]
        if (!is.null(x.names))
           x.names <- x.names[opt.order]
    }
    if (scale.percent) {
        colM <- apply(d, 2, max)
        colM <- floor((colM + 5)/5) * 5
        colM[colM < min.width] <- min.width
        colM.sum <- sum(colM)
    }
    else {
        colM.sum <- sum(graph.widths)
        colM <- graph.widths
    }
    xLen <- xRight - xLeft
    xInc <- xLen - ((nsp + 1) * xSpace)
    inc <- xInc/colM.sum
    if (inc < 0.0)
       stop("Too many variables, curves will be too small. Reduce xSpace.")
    x1 <- xLeft
#    par(fig = c(x1, x1+0.4, yStart, yTop))
    par(fig = figCnvt(orig.fig, c(x1, min(x1+0.4, .9), yBottom, yTop)), new=add)
    if (y.rev) {
        tmp <- ylim[1]
        ylim[1] <- ylim[2]
        ylim[2] <- tmp
    }
    plot(0, 0, cex = 0.5, xlim = c(0, 1), axes = FALSE, type = "n", yaxs = "r", ylim = ylim, ...)
    usr1 <- par("usr")
    if (y.axis) {
      if (is.null(y.tks))
         y.tks <- axTicks(2)
      ax <- axis(side = 2, las = 1, at = y.tks, labels = as.character(y.tks), cex.axis=cex.yaxis, xpd=NA)
      x1 <- x1 + xSpace
      mtext(title, adj = 0, line = 5, cex = cex.title)
      mtext(ylabel, side = 2, line = 2.5, cex=cex.ylabel)
    }
    ty <- ifelse(plot.line, "l", "n")
    tcll <- -.3
    if ("tcl" %in% names(fcall))
       tcll <- eval(fcall$tcl)
    spc <- 0
    if ("las" %in% names(fcall)) {
       if ((eval(fcall$las)) == 2)
          spc = 0.3
    }
    figs <- vector("list", length=nsp)
    usrs <- vector("list", length=nsp)
    for (i in 1:nsp) {
        par(new = TRUE)
        par(lend = "butt")
        if (scale.percent) {
            inc2 <- inc * colM[i]
            par(fig = figCnvt(orig.fig, c(x1, x1 + inc2, yBottom, yTop)))
            plot(0, 0, cex = 0.5, xlim = c(0, colM[i]), axes = FALSE, 
                xaxs = "i", type = "n", yaxs = "r", ylim = ylim, xlab="", ylab="", ...)
            if (plot.symb) {
               points(d[, i], yvar, pch=symb.pch, cex=symb.cex, xpd=NA)
            }
            if (plot.poly) {
               y <- c(min(yvar, na.rm=TRUE), yvar, max(yvar, na.rm=TRUE))
               x <- c(0, d[, i], 0)
               if (exag[i]) {
                 x2 <- c(0, d[, i]*exag.mult[i], 0)
                 polygon(x2, y, col = col.exag[i], border = NA)
               }
               polygon(x, y, col = cc.poly[i], border = cc.poly.line[i], lwd=lwd.poly)
            }
            if (is.logical(plot.bar)) {
               if (plot.bar) {
                  if (sep.bar) {
                      segments(rep(0, nsam), yvar, d[, i], yvar, lwd = lwd.bar, col = cc.bar)
                  } else {
                   segments(rep(0, nsam), yvar, d[, i], yvar, lwd = lwd.bar, col = cc.bar[i])
                  }
               }
            } else {
              if (plot.bar=="Full") {
                abline(h=yvar, col=cc.bar)
              }
            }
            lines(c(0, 0), c(min(yvar, na.rm=TRUE), max(yvar, na.rm=TRUE)), ...)
            if (ty == "l") 
               lines(d[, i], yvar, col = cc.line[i], lwd = lwd.line)
            if (add.smooth[i]) {
              tmp <- data.frame(x=yvar, y=d[, i])
              tmp <- na.omit(tmp)
              lo <- lowess(tmp, f=smooth.span[i])
              lines(lo$y, lo$x, col=col.smooth, lwd=lwd.smooth)
            }
            xlabb <- seq(0, colM[i], by = x.pc.inc[i])
            if (x.pc.lab) {
               xlabbt <- as.character(xlabb)
               if (x.pc.omit0)
                  xlabbt[1] <- ""
               mgpX <- if (is.null(mgp)) { c(3,max(0.0, spc-tcll),0) } else { mgp }
               axis(side = 1, at = xlabb, labels = xlabbt, mgp=mgpX, cex.axis=cex.axis, ...)
            }
            else
               axis(side = 1, at = xlabb, labels = FALSE, ...)
            x1 <- x1 + inc2 + xSpace
        }
        else {
            inc2 <- inc * colM[i]
            par(fig = figCnvt(orig.fig, c(x1, min(1, x1 + inc2), yBottom, yTop)))
            if (!is.null(minmax)) {
               plot(d[, i], yvar, cex = 0.5, axes = FALSE, xaxs = "i", 
                type = "n", yaxs = "r", ylim = ylim, xlim=c(minmax[i, 1], minmax[i,2]), ...)
            } else {
               plot(d[, i], yvar, cex = 0.5, axes = FALSE, xaxs = "i", 
                type = "n", yaxs = "r", ylim = ylim, ...)
            }
            tks <- axTicks(1)
            us <- par("usr")
            if (plot.symb) {
               points(d[, i], yvar, pch=symb.pch, cex=symb.cex, xpd=NA)
            }
            if (plot.poly) {
               y <- c(min(yvar, na.rm=TRUE), yvar, max(yvar, na.rm=TRUE))
               x <- c(us[1], d[, i], us[1])
               if (exag[i]) {
                 x2 <- c(0, d[, i]*exag.mult[i], 0)
                 polygon(x2, y, col = col.exag[i], border = NA)
               }
               polygon(x, y, col = cc.poly[i], border = cc.poly.line[i], lwd=lwd.poly)
            }
            if (is.logical(plot.bar)) {
              if (plot.bar) {
               if (sep.bar) {
                 segments(rep(us[1], nsam), yvar, d[, i], yvar, lwd = lwd.bar, col = cc.bar)
               } else {
                 segments(rep(us[1], nsam), yvar, d[, i], yvar, lwd = lwd.bar, col = cc.bar[i])               }
               }
            } else {
               if (plot.bar=="Full") {
                  abline(h=yvar, col=cc.bar)
               }
            }
            lines(c(us[1], us[1]), c(min(yvar, na.rm=TRUE), max(yvar, na.rm=TRUE)), 
                ...)
            if (ty == "l") 
                lines(d[, i], yvar, col = cc.line[i], lwd = lwd.line)
            if (add.smooth[i]) {
              tmp <- data.frame(x=yvar, y=d[, i])
              tmp <- na.omit(tmp)
              lo <- lowess(tmp, f=smooth.span[i])
              lines(lo$y, lo$x, col=col.smooth, lwd=lwd.smooth)
            }
            mgpX <- if (is.null(mgp)) { c(3, max(0.0, spc-tcll),0) } else { mgp }
            if (scale.minmax) {
               nn <- length(axTicks(1))
               tk <- c(axTicks(1)[1], axTicks(1)[nn])
               axis(side = 1, at = tk, labels = as.character(tk), cex.axis=cex.axis, mgp=mgpX, ...)
            }
            else {
               axis(side = 1, cex.axis=cex.axis, mgp=mgpX, ...)
            }
            x1 <- x1 + inc2 + xSpace
        }
#        tks1 <- axTicks(1)
        usr2 <- par("usr")
        tks1 <- usr2[1]

        r <- (usr1[4] - usr1[3]) * 0.01
        pos <- usr1[4]+r
        if (y.rev)
           pos <- usr1[4]-r
        if (srt.xlabel < 90)
           text(tks1[1], pos, labels=x.names[i], adj = c(0, 0), srt=srt.xlabel, cex = cex.xlabel, xpd=NA)
        else
           text(tks1[1], pos, labels=x.names[i], adj = c(0, 1), srt=srt.xlabel, cex = cex.xlabel, xpd=NA)
           
        usrs[[i]] <- usr2   
        figs[[i]] <- par("fig")
    }
    if (!is.null(clust)) {
        par(fig = figCnvt(orig.fig, c(x1, xRight+clust.width, yBottom, yTop)))
        par(mar=c(0,0,0,0))
        par(new = TRUE)
#        plot(clust, horiz = TRUE, xaxt.rev=yaxt.rev, leaflab = "none", cex.axis = 0.5, yaxt.rev=TRUE)
#        if(y.rev)
#           clust <- rev(clust)
        mgpX <- if (is.null(mgp)) { c(2, .5, 0) } else { mgp }
        plot(clust, xvar=yvar, horiz=TRUE, x.rev=y.rev, labels=rep("", length(yvar)), hang=-1, mgp=mgpX, cex.axis=cex.axis, ...)
    }
    par(mai = oldmai)
    oldfig[oldfig < 0] <- 0
    par(fig = oldfig)
    ll <- list(call=fcall, box=c(xLeft=xLeft, xRight=xRight, yBottom=yBottom, yTop=yTop), usr = usr1, yvar=yvar, ylim=ylim, y.rev=y.rev, figs=figs, usrs=usrs)
    invisible(ll)
}

addZone <- function(x, upper, lower=NULL, ...) {
   oldpar <- par(c("fig", "mar", "usr"))
   par(fig=x$box)
   par(mar=c(0,0,0,0))
   par(usr=c(0, 1, x$usr[3], x$usr[4]))
   if (is.null(lower))
      segments(0, upper, 1, upper, xpd=NA, ...)
   else
      rect(0, lower, 1, upper, ...)
   par(oldpar)
}

addClustZone <- function(x, clust, nZone, ...) {
   oldpar <- par(c("fig", "mar", "usr"))
   par(fig=x$box)
   par(mar=c(0,0,0,0))
   par(usr=c(0, 1, x$usr[3], x$usr[4]))
   cc <- cutree(clust, k=nZone)
   zn <- which(diff(cc)>0)
#   if (x$yaxt.rev)
#      x$yvar <- rev(x$yvar)
   zone <- (x$yvar[zn] + x$yvar[zn+1]) / 2
   segments(0, zone, 1, zone, xpd=NA, ...)
   par(oldpar)
}

strat.plot.simple <- function(y1, x1, y2=NULL, x2=NULL, col=c("blue", "red"), sort.vars=c("original","wa", "alphabetical"), ylim=range(x1), y.rev=FALSE, type=c("b", "l"), subset=c(1:ncol(y1)), ...) {
   nsp.y1 <- ncol(y1)
   y1 <- y1[, subset]
   if (ncol(y1) > 50)
      stop("You have more than 50 columns in the species data, split the data into smaller subsets")
   sort.vars <- match.arg(sort.vars)
   if (sort.vars == "original") {
      or1 <- order(colnames(y1))
      ord <- order(or1)
   } else {
      if (sort.vars == "wa") {
          or1 <- order(colnames(y1))
          wa.sc <- apply(y1[, or1], 2, function(x, env) { sum(x*env, na.rm=TRUE) / sum(x, na.rm=TRUE) }, env=x1)
          ord <- order(wa.sc)
      }
      else {
         ord <- 1:ncol(y1)
      }
   }
#   require(lattice)
   s <- stack(y1)
   s$x <- rep(x1, times=ncol(y1))
   s$set <- 1
   if (!is.null(y2) | !is.null(x2)) {
     if (nsp.y1 != ncol(y2))
        stop("Number of columns different in y1 and y2")
     y2 <- y2[, subset] 
     s2 <- stack(data.frame(y2))
     s2$x <- rep(x2, times=ncol(y2))
     s2$set <- 2
     s <- rbind(s, s2)
   }
   if (y.rev)
      ylim <- rev(ylim)
   xyplot(x ~ values| ind, data = s, type=type, distribute.type=TRUE, groups = s$set, col=c("blue", "red"), ylim=ylim, index.cond=list(ord), ylab="", xlab="", ...)
}
