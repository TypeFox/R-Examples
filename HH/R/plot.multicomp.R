plot.multicomp.hh <-
  function(x, ylabel = x$ylabel, href = 0, uniform = TRUE,
           plt.in = c(0.2, 0.9, 0.1, 0.9),
           ##         fig.in=par("fig"),
           x.label.adj=1,
           xrange.include=href,
           xlim,
           comparisons.per.page=21,
           col.signif=1, col.not.signif=1,
           lty.signif=4, lty.not.signif=4,
           lwd.signif=1, lwd.not.signif=1,
           ...,
           xlabel.print=TRUE, y.axis.side=2, ylabel.inside=FALSE)
{
  if.R(r=
       NextMethod("plot")
       ,
       s={
         table <- x$table
         method <- x$method
         if(method == "lsd")
           method <- "LSD"
         if(method == "scheffe")
           method <- "Scheffe"
         if(method == "bon")
           method <- "Bonferroni"
         if(method == "sidak")
           method <- "Sidak"
         if(method == "tukey")
           method <- "Tukey"
         if(method == "dunnett")
           method <- "Dunnett"
         if(method == "sim")
           method <- "simulation-based"
         if(method == "tpmc")
           method <- "Cheung and Chan"
         if(x$error.type == "fwe")
           xlabel <- paste("simultaneous ", as.character(100 * (1 - x$
                                                                alpha)), "% confidence limits,", method, "method")
         if(x$error.type == "cwe")
           xlabel <- paste("individual", as.character(100 * (1 - x$alpha)),
                           "% confidence limits,", method, "method")

         old.par <- par()       ## this is imperative.
         ##  There is a bug in plt through at least S-Plus 7.0.6
         on.exit({par(new.usr); par(new=FALSE)})  ## this is imperative.
         old.plt <- par(plt = plt.in)
         ##        par(plt = plt.in, fig=fig.in)
         ##        old.par <- par()      ## at least one other par value is affected
         ##        par(new=FALSE)        ## new seems to be reset when fig is set
         ##        on.exit({par(old.par)
         ##                 par(new=FALSE)
         ##                 par(new.usr)}) ## restore everything except new and usr

         ptot <- nrow(table)
         cpp <- comparisons.per.page
         npage <- ceiling(ptot/cpp)
         for(page in 1:npage) {
           estimate <- table[((page - 1) * cpp + 1):min(((page - 1) * cpp +
                                                         cpp), ptot), "estimate"]
           lower <- table[((page - 1) * cpp + 1):min(((page - 1) * cpp + cpp),
                                                     ptot), "lower"]
           upper <- table[((page - 1) * cpp + 1):min(((page - 1) * cpp + cpp),
                                                     ptot), "upper"]
           labels <- dimnames(table)[[1]][((page - 1) * cpp + 1):min(((
                                                                       page - 1) * cpp + cpp), ptot)]
           signif <- (lower * upper > 0) ## hh
           p <- length(estimate)
           lower <- ifelse(lower == ( - Inf), NA, lower)
           upper <- ifelse(upper == Inf, NA, upper)
           if (missing(xlim)) {
             if(page == 1 || !uniform) {
               xmin <- min(c(lower, upper, estimate, href,
                             xrange.include), na.rm = TRUE)
               xmax <- max(c(lower, upper, estimate, href,
                             xrange.include), na.rm = TRUE)
               xrange <- range(pretty(c(xmin, xmax), 10))
             }}
           else {
             xmin <- xlim[1]
             xmax <- xlim[2]
             xrange <- range(pretty(c(xmin, xmax), 10))
           }
           lines.lower <- ifelse(is.na(lower), min(xrange), lower)
           lines.upper <- ifelse(is.na(upper), max(xrange), upper)
           ends.lower <- !is.na(lower)
           ends.upper <- !is.na(upper)
           heights <- cpp:(cpp - p + 1)
           plot(estimate, heights, ylab = "", ylim = c(0, cpp+1), yaxt = "n",
                xlim = xrange, xlab = "", xaxt = "n", pch = 16, bty =
                "n", ...,
                type="n")
           new.usr <- par()["usr"]
           if (sum(signif))
             points(estimate[signif], heights[signif], pch=16, col=col.signif)
           if (sum(!signif))
             points(estimate[!signif], heights[!signif], pch=16, col=col.not.signif)

           xlower <- lower[ends.lower & signif]
           if(length(xlower) > 0)
             points(xlower, heights[ends.lower & signif], pch = "(",
                    col=col.signif)
           xupper <- upper[ends.upper & signif]
           if(length(xupper) > 0)
             points(xupper, heights[ends.upper & signif], pch = ")",
                    col=col.signif)
           segments(lines.lower[signif], heights[signif],
                    lines.upper[signif], heights[signif],
                    lty = lty.signif, lwd = lwd.signif,
                    col=col.signif)

           xlower <- lower[ends.lower & !signif]
           if(length(xlower) > 0)
             points(xlower, heights[ends.lower & !signif], pch = "(",
                    col=col.not.signif)
           xupper <- upper[ends.upper & !signif]
           if(length(xupper) > 0)
             points(xupper, heights[ends.upper & !signif], pch = ")",
                    col=col.not.signif)
           segments(lines.lower[!signif], heights[!signif],
                    lines.upper[!signif], heights[!signif],
                    lty = lty.not.signif, lwd = lwd.not.signif,
                    col=col.not.signif)

           if(!is.null(href))
             segments(href, rep(cpp+1 - p - 1, length(href)), href,
                      rep(cpp+1, length(href)))
           if (sum(signif)) {
             ##text(xmax + 0.1 * (xmax - xmin),
             if (ylabel.inside) {
               ## axis(side=y.axis.side,
               ##     at=heights[signif], labels=FALSE, col=col.signif,
               ##     adj = x.label.adj, ticks=FALSE, col.axis=col.signif)
               if (y.axis.side==2) {
                 axis(side=4, pos=xrange[1],
                      at=heights[signif], labels=labels[signif], col=col.signif,
                      adj = x.label.adj, ticks=FALSE, col.axis=col.signif)
               }
               if (y.axis.side==4) {
                 axis(side=2, pos=xrange[2]-diff(xrange)/25,
                      at=heights[signif], labels=labels[signif], col=col.signif,
                      adj = x.label.adj, ticks=FALSE, col.axis=col.signif)
               }
             }
             else
               axis(side=y.axis.side,
                    at=heights[signif], labels=labels[signif], col=col.signif,
                    adj = x.label.adj, ticks=FALSE, col.axis=col.signif)
           }
           if (sum(!signif)) {
             ##text(xmax + 0.1 * (xmax - xmin),
             if (ylabel.inside) {
               if (y.axis.side==2) {
                 axis(side=4, pos=xrange[1],
                      at=heights[!signif], labels=labels[!signif], col=col.not.signif,
                      adj = x.label.adj, ticks=FALSE, col.axis=col.not.signif)
               }
               if (y.axis.side==4) {
                 axis(side=2, pos=xrange[2]-diff(xrange)/25,
                      at=heights[!signif], labels=labels[!signif], col=col.not.signif,
                      adj = x.label.adj, ticks=FALSE, col.axis=col.not.signif)
               }
             }
             else
               axis(side=y.axis.side,
                    at=heights[!signif], labels=labels[!signif], col=col.not.signif,
                    adj = x.label.adj, ticks=FALSE, col.axis=col.not.signif)
           }
       axis(1, at = pretty(c(xmin, xmax), 10), pos = cpp - p)
           if (xlabel.print) {
             text((xmax + xmin)/2, (cpp - p - 3)*(22/cpp)-(cpp/22),
                  xlabel, adj = 0.5)
             text((xmax + xmin)/2, (cpp - p - 3)*(22/cpp)-2*(cpp/22),
                  paste("response variable:", ylabel), adj = 0.5)
           }
           lines(c(xrange[1], xrange[1]), c(cpp - p, cpp+1))
           lines(c(xrange[2], xrange[2]), c(cpp - p, cpp+1))
           lines(c(xrange[1], xrange[2]), c(cpp+1, cpp+1))
           lines(c(xrange[1], xrange[2]), c(cpp - p, cpp - p))
         }
       }
       )
}


plot.matchMMC <- function(...)
  .Defunct("plotMatchMMC", package="HH")

plotMatchMMC <- function(x, ...,
                         xlabel.print=FALSE,
                         cex.axis=par()$cex.axis,
                         col.signif='red', main="",
                         ylabel.inside=FALSE,
                         y.axis.side=4,
                         adjusted=FALSE) {
  if (inherits(x, "mmc.multicomp")) {
    x <- if (!is.null(x$lmat)) x$lmat else x$mca
  }
  if (!inherits(x, "multicomp.hh"))
    stop('plotMatchMMC requires a "multicomp.hh" or "mmc.multicomp" object.')
  if.R(s={
    old.xpd <- par(xpd=TRUE)
    xlim <- par()$usr[1:2]
    xlim <- xlim + c(.1,-.1)*diff(xlim)
    plot(x, xrange.include=xlim, xaxs="d", ## xaxs="i",
         main=main, xlab="", ...,
         col.signif=col.signif, lty.signif=1, xlabel.print=xlabel.print,
         plt=par()$plt+c(0,0,-.25,.05), x.label.adj=0,
         y.axis.side=y.axis.side, ylabel.inside=ylabel.inside)
    invisible(par(old.xpd))
  },
       r={
         if (adjusted)
           plot.multicomp.adjusted(x,
                xlim=par()$usr[1:2], xaxs="i", yaxt="n",
                main=main, xlab="", cex.axis=cex.axis)
         else
           plot(x,
                xlim=par()$usr[1:2], xaxs="i", yaxt="n",
                main=main, xlab="", cex.axis=cex.axis)
         signif <- apply(x$table[,c("lower","upper"), drop=FALSE], 1, prod) > 0
         yval <- rev(seq(along=signif))
         if (!all(!signif)) {
           if (ylabel.inside)
             mtext(names(signif)[signif], at=yval[signif],
                   side=4, adj=1, line=-.6, col=col.signif, las=1)
           else
             axis(4, at=yval[signif], labels=names(signif)[signif],
                  col=col.signif, col.axis=col.signif,
                  las=1, tck=-.01, mgp=c(3,.5,0), cex.axis=cex.axis)
           lower <- x$table[signif, "lower", drop=FALSE]
           upper <- x$table[signif, "upper", drop=FALSE]
           lower[lower==-Inf] <- par()$usr[1]
           upper[upper==Inf] <- par()$usr[2]
           segments(lower, yval[signif],
                    upper, yval[signif],
                    col=col.signif)
         }
         if (!all(signif)) {
           if (ylabel.inside)
             mtext(names(signif)[!signif], at=yval[!signif],
                   side=4, adj=1, line=-.6, las=1)
           else
             axis(4, at=yval[!signif], labels=names(signif)[!signif],
                  las=1, tck=-.01, mgp=c(3,.5,0), cex.axis=cex.axis)
         }
       }
       )
}
## assignInNamespace("plotMatchMMC", plotMatchMMC, "HH")
## source("c:/HOME/rmh/HH-R.package/HH/R/plot.multicomp.R")
