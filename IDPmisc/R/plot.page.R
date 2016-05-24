## plot.page.R
## helper function for longtsPlot

plot.page <- function(y1, y2, names1, names2,
                      startP, upf, fpp, overlap,
                      x.at, x.ann, x.tick,
                      y1.at, y1.ann, y1.tick,
                      y2.at, y2.ann, y2.tick,
                      ny.ann, cex.ann,
                      xlab, y1lab, y2lab, las, col.y1, col.y2, cex.lab,
                      y1lim, y2lim,
                      lty1, lty2, lwd1, lwd2, col1, col2,
                      leg, y1nam.leg, y2nam.leg,
                      ncol.leg, cex.leg=1.5,
                      h1, h2, col.h1, col.h2,
                      mgp, main, cex.main, xpd, cex,
                      type1, type2, pch1, pch2, cex.pt1, cex.pt2){
    ## internal function
    ## Author:  Rene Locher
    ## Version: 2011-10-07

    on.exit(options(warn=options()$warn, add=TRUE))
    options(warn=-1)

    ## preparing arguments for plotting depending on type
    Nmax <- max(length(type1),length(pch1),
                length(lwd1),length(lty1),
                length(cex.pt1),length(col1))
    type1 <- rep(type1,Nmax)[1:Nmax]
    lwd1 <- rep(lwd1,Nmax)[1:Nmax]
    lty1 <- rep(lty1,Nmax)[1:Nmax]
    pch1 <- rep(pch1,Nmax)[1:Nmax]
    cex.pt1 <- rep(cex.pt1,Nmax)[1:Nmax]
    col1 <- rep(col1,Nmax)[1:Nmax]

    Nmax <- max(length(type2),length(pch2),
                length(lwd2),length(lty2),
                length(cex.pt2),length(col2))
    type2 <- rep(type2,Nmax)[1:Nmax]
    lwd2 <- rep(lwd2,Nmax)[1:Nmax]
    lty2 <- rep(lty2,Nmax)[1:Nmax]
    pch2 <- rep(pch2,Nmax)[1:Nmax]
    cex.pt2 <- rep(cex.pt2,Nmax)[1:Nmax]
    col2 <- rep(col2,Nmax)[1:Nmax]

    for (ff in 0:(fpp-1)) {
        st <- startP+ff*upf
        ee <- startP+(ff+1)*upf+overlap
        st1 <- start(y1)
        st1 <- max(st,st1[1]+(st1[2]-1)/frequency(y1))
        ee1 <- end(y1)
        ee1 <- min(ee,ee1[1]+(ee1[2]-1)/frequency(y1))

        if (st1 < ee1) {
            opar <- par(xpd=xpd)
            err <- try(plot(window(y1, start = st1, end = ee1),
                            col = col1, lty = lty1, lwd = lwd1,
                            type = "n", pch = pch1, cex = cex,
                            ylim = y1lim, xlim=c(st, ee),
                            plot.type = "single", ## see plot.ts!!
                            xlab = xlab, ylab="", axes = FALSE),
                       silent=TRUE)
            par(opar)
            if (!is.null(err)) {
                if (regexpr("margins too large",geterrmessage())>0) {
                    stop("\nToo many figures per page defined.\nChoose smaller value for 'fpp' and use slide=TRUE\n", call. = FALSE)} else
                stop(paste("\nUndefined error:\n",
                           geterrmessage(),sep=""), call. = FALSE)
            }

            for (ii in 1:ncol(y1)) {
                lines(window(y1[,ii], start = st1, end = ee1),
                      col = col1[ii], lty = lty1[ii], lwd = lwd1[ii],
                      pch = pch1[ii], cex = cex.pt1[ii],
                      type=type1[ii])
            }

            if (!is.null(h1)) abline(h=h1, col=col.h1, xpd=FALSE)

            if (!is.null(y2)) {
                st2 <- start(y2)
                st2 <- max(st,st2[1]+(st2[2]-1)/frequency(y2))
                ee2 <- end(y2)
                ee2 <- min(ee,ee2[1]+(ee2[2]-1)/frequency(y2))

                if (st2<ee2) {
                    par(xpd=xpd)
                    for (ii in 1:ncol(y2)) {
                        lines(window(y2[,ii],start=st2,end=ee2),
                              col=col2[ii], lty=lty2[ii], lwd=lwd2[ii],
                              type=type2[ii], pch=pch2[ii],
                              cex=cex.pt2[ii])
                    }
                    par(opar)
                }
                if (!is.null(h2)) abline(h=h2, col=col.h2, xpd=FALSE)
            }

            if (!is.null(x.tick))
                axis(1,at=x.tick,labels=FALSE,tcl=-0.3, xpd=FALSE, las=las)
            axis(1, at = x.at, labels = x.ann,
                 cex.axis = cex.ann, tcl = -0.5, xpd = FALSE, las=las)

            if (!is.null(y1.tick))
                axis(2, y1.tick, labels = FALSE, col.axis = col.y1,
                     tcl=-0.3, xpd=FALSE, las=las)
            axis(2, at = y1.at, labels = y1.ann, col.lab = col.y1,
                 col.axis = col.y1, cex.axis = cex.ann, tcl=-0.5,
                 xpd = FALSE, las=las)
            mtext(text = y1lab, side = 2, line = mgp[1],
                  col = col.y1, cex = cex.lab)

            if (!is.null(y2)) {
                if (!is.null(y2.tick))
                    axis(4, y2.tick, labels=FALSE, col.axis=col.y2,
                         tcl=-0.3, xpd=FALSE, las=las)
                axis(4, at=y2.at, labels=y2.ann, col.lab=col.y2,
                     col.axis=col.y2, cex.axis=cex.ann, tcl=-0.5, xpd=FALSE,
                     las=las)
                mtext(text=y2lab, side=4, line=mgp[1], col=col.y2,
                      cex=cex.lab)
            }
            box()
        } else frame()
    }

    if (leg) {

        ## prepatation of lwd, lty and pch for legend
        iNA <- type1 == "p"
        lwd1[iNA] <- lty1[iNA] <- NA

        iNA <- type2 == "p"
        lwd2[iNA] <- lty2[iNA] <- NA

        pch1[type1 %in% c("l","s","S")] <- NA
        pch2[type2 %in% c("l","s","S")] <- NA

        par(xpd=NA)
        plot(0:1, 0:1, type="n", an = FALSE, axes=FALSE)
        cH <- strheight("A")
        cW <- strwidth("A")

        if (is.null(y2)) { ## right axis does not exist
            if (is.null(ncol.leg)) ncol.leg <- ncol(y1)
            le <- legend(0.5,0.5, legend=names1,
                         col=col1, lwd=lwd1, lty=lty1,
                         pch=pch1, pt.cex=cex.pt1,
                         ncol=ncol.leg, xjust=0.5, yjust=0.5, cex=cex.leg,
                         bty="n",plot=FALSE)
            if (is.null(y1nam.leg)) ## without axis title in legend
                legend(0.5, le$rect$h,
                       legend=names1,
                       col=col1, lwd=lwd1, lty=lty1,
                       ncol=ncol.leg,
                       xjust=0.5, yjust=1, cex=cex.leg,
                       bty="n") else { ## with axis title in legend
                           le <- legend(0.5, le$rect$h,
                                        legend=names1,
                                        col=col1, lwd=lwd1, lty=lty1,
                                        pch=pch1, pt.cex=cex.pt1,
                                        ncol=ncol.leg,
                                        xjust=0.5, yjust=1, cex=cex.leg,
                                        bty="n")
                           text(0.5,le$rect$h+cH*cex.leg,
                                y1nam.leg, cex=cex.leg, adj=0.5)
                       }
        } else { ## right axis exists
            if (is.null(ncol.leg)) ncol.leg <- max(ncol(y1),ncol(y2))
            if (is.null(y1nam.leg)) y1nam.leg <- "left axis:"
            if (is.null(y2nam.leg)) y2nam.leg <- "right axis:"
            text.w <- max(sapply(c(names1,names2),strwidth, cex=cex.leg))

            le1 <- legend(0.5,0.5,legend=names1,
                          col=col1,lwd=lwd1,lty=lty1,
                          pch=pch1, pt.cex=cex.pt1,
                          ncol=ncol.leg, xjust=0.5, yjust=1, cex=cex.leg,
                          bty="n", text.width=text.w, plot=FALSE)
            le2 <- legend(0.5,0.5, legend=names2,
                          col=col2,lwd=lwd2,lty=lty2,
                          pch=pch2, pt.cex=cex.pt2,
                          ncol=ncol.leg, xjust=0.5, yjust=1, cex=cex.leg,
                          bty="n", text.width=text.w, plot=FALSE)

            leg.w <- min(max(le1$rect$w,le2$rect$w),1)

            le1 <- legend(0.5-leg.w/2, le1$rect$h+le2$rect$h-2*cH,
                          legend=names1,
                          col=col1, lwd=lwd1, lty=lty1,
                          pch=pch1, pt.cex=cex.pt1,
                          ncol=ncol.leg, xjust=0, yjust=1, cex=cex.leg,
                          bty="n", text.width=text.w)
            le2 <- legend(0.5-leg.w/2, le2$rect$h-cH,
                          legend=names2,
                          col=col2, lwd=lwd2, lty=lty2,
                          pch=pch2, pt.cex=cex.pt2,
                          ncol=ncol.leg, xjust=0, yjust=1, cex=cex.leg,
                          bty="n", text.width=text.w)

            text(0.5-(leg.w+cW*cex.leg)/2, le1$text$y[1], y1nam.leg,
                 cex=cex.leg, adj=1, col=col.y1)
            text(0.5-(leg.w+cW*cex.leg)/2, le2$text$y[1], y2nam.leg,
                 cex=cex.leg, adj=1, col=col.y2)

        }
        par(xpd=xpd)
    }
    if (!is.null(main)) mtext(text=main,side=3,line=0,outer=T,cex=cex.main)
} #plot.page
