plot.mmc.multicomp <-
  function(x,
           xlab="contrast value",
           ylab=none$ylabel,
           focus=none$focus,
           main= main.method.phrase,
           main2=main2.method.phrase,
           main.method.phrase= paste("multiple comparisons of means of", ylab),
           main2.method.phrase=paste("simultaneous ",
             100*(1-none$alpha),"% confidence limits, ",
             method, " method", sep="" ),
           ry.mmc=TRUE,
           key.x=par()$usr[1]+ diff(par()$usr[1:2])/20,
           key.y=par()$usr[3]+ diff(par()$usr[3:4])/3,
           method=if (is.null(mca)) lmat$method else mca$method,
           print.lmat=(!is.null(lmat)),
           print.mca=(!is.null(mca) && (!print.lmat)),
           iso.name=TRUE,
           x.offset=0,
           col.mca.signif='red',  col.mca.not.signif='black',
           lty.mca.signif=1,  lty.mca.not.signif=6,
           lwd.mca.signif=1,  lwd.mca.not.signif=1,
           col.lmat.signif='blue', col.lmat.not.signif='black',
           lty.lmat.signif=1, lty.lmat.not.signif=6,
           lwd.lmat.signif=1, lwd.lmat.not.signif=1,
           lty.iso=7, col.iso='darkgray', lwd.iso=1,
           lty.contr0=2, col.contr0='darkgray', lwd.contr0=1,
           decdigits.ybar=2,
           ...
           ) {
    if (inherits(x,"mmc.multicomp")) {
      lmat <- x$lmat
      none <- x$none
      mca  <- x$mca
    }
    else
      stop("x must be an mmc.multicomp object")
    if (is.null(mca) && is.null(lmat))
      stop("at least one of mca and lmat are required in x")
    if (is.null(none))
      stop("none is required in x")


    ## if (missing(main)&& !is.null(x$main)) main  <- x$main ## easy to read
    where.is.main <- match("main", names(x), 0)
    if (missing(main)  && where.is.main)  main  <- x$main ## work correctly
    if (missing(main2) && !is.null(x$main2)) main2 <- x$main2
    
    ## set up plot, using ybar from none
    ybar <- none$table[,"estimate"]

    ry <- range(ybar) + c(-.05,0)*diff(range(ybar))
    if (is.numeric(ry.mmc)) ry <- ry.mmc
    else if (ry.mmc && !is.null(x$ry)) ry <- x$ry

    old.par <- if.R(r=par(no.readonly=TRUE), s=par()) ## save par() settings
   
    if (par()$xaxs=="d" || par()$yaxs=="d")
      usr.s <- par()$usr
    else
      { ## force a square plotting region rotated 45 degrees
        if.R(r={pin.m <- par()$pin
                plt.m <- par()$plt
                usr.m <- par()$usr},
             s=par(pty="s"))
        plot(ylim=ry,
             xlim=c(-1,1) * diff(ry) + x.offset,
             x=c(0,0), y=ry,
             type="n", xlab="", ylab="", axes=FALSE,
             asp=2)   ## aspect honored by R, ignored by S-Plus
        
        usr.s <- par()$usr
        pin.s <- par()$pin ## width and height of plot, measured in inches
        plt.s <- par()$plt ## coordinates of the plot region as fractions
                           ## of the current figure region
        
        if.R(s={par(pty="m")
                pin.m <- par()$pin
                usr.m.new <- c(usr.s[1:2] * (pin.m/pin.s)[1], usr.s[3:4])
                par(usr=usr.m.new)
                par(new=TRUE)},
             r={})
##        if.R(r={par(plt=plt.m,
##                    pin=pin.m)},
##             s={})
      }

    if.R(s=
         plot(y=ry,
              x=c(0,0),
              type="n", xaxs="d", yaxs="d",
              yaxt="n",
              xlab=xlab, ylab="", main=main)
         ,r={
##          plot(y=ry,
##               x=c(0,0),
##               type="n", ## xaxs="d", yaxs="d",  ## axis style "d" unimplemented in R
##               xlim=par()$usr[1:2], ylim=par()$usr[3:4],  ## axis style "d" unimplemented in R
##               yaxt="n",
##               xlab=xlab, ylab="", main=main)
           box()
           title(main=main, xlab=xlab)
           axis(1)
         }
         )    
    mtext(side=3, main2, line=.5, cex=1.2)
    mtext(side=2, "mean", line=4, adj=0,
          at=par()$usr[3]-.01*diff(par()$usr[3:4]), las=1)
    mtext(side=2, ylab, line=4, adj=0,
          at=par()$usr[3]-.05*diff(par()$usr[3:4]), las=1)
    mtext(side=2, paste(focus, "level"), line=2, adj=0,
          at=par()$usr[3]-.11*diff(par()$usr[3:4]), las=1)
    mtext(side=4, "contrast", line=0, adj=.5,
          at=par()$usr[3]-.11*diff(par()$usr[3:4]), las=1)

    ## the additional feature of the Hsu Peruggia plot is the
    ## iso-lines for each level of the factor
    yb <- as.vector(ybar)
    names(yb) <- seq(along=yb)
    ybo <- order(yb)
    yb <- yb[ybo]
    yb.min <- min(yb)
    yb.max <- max(yb)
    segments(yb.min-yb, (yb.min+yb)/2,
             yb.max-yb, (yb.max+yb)/2, lty=lty.iso, col=col.iso, lwd=lwd.iso)
    segments(yb-yb.min, (yb+yb.min)/2,
             yb-yb.max, (yb+yb.max)/2, lty=lty.iso, col=col.iso, lwd=lwd.iso)
    abline(v=0, lty=lty.contr0, col=col.contr0, lwd=lwd.contr0)

    axis(1, at=yb.min-yb, tck=.02, labels=FALSE)
    ##  axis(1, at=yb.min-yb, tick=FALSE, labels=names(yb), line=-2.5)
    axis(1, at=yb-yb.min, tck=.02, labels=FALSE)
    ##  axis(1, at=yb-yb.min, tick=FALSE, labels=names(yb), line=-2.5)

    yb.names <- if (iso.name) names(ybar)[ybo] else names(yb)
    yb.min.index <- match(yb.min,yb)
    text(yb.min-yb[-yb.min.index], (yb.min+yb)[-yb.min.index]/2-.025*diff(usr.s[3:4]), yb.names[-yb.min.index], adj=.9, col=col.iso)
    text(yb[-yb.min.index]-yb.min, (yb+yb.min)[-yb.min.index]/2-.025*diff(usr.s[3:4]), yb.names[-yb.min.index], adj=.1, col=col.iso)
    text(0,                        (yb+yb.min)[ yb.min.index]/2-.025*diff(usr.s[3:4]), yb.names[ yb.min.index], adj=.5, col=col.iso)

    ## y-axis labels
    for (i in seq(along=ybar)) {
      axis(2, at=ybar[i], adj=1,                   ## exterior labels and ticks
           labels=format(round(ybar[i], decdigits.ybar), nsmall=decdigits.ybar),
           las=1)
      ##   both format and round are needed!
      left.interior.label.args <-
        list(2, at=ybar[i], tick=FALSE,
             line=-2, adj=0, las=1,
             label=
             if (iso.name) {
               names(ybar[i])}
             else {
               paste(i, ": ", names(ybar[i]), sep="")})
      if.R(r=names(left.interior.label.args)[5] <- "hadj",
           s={})
      do.call("axis", left.interior.label.args)
      axis(2, at=ybar[i], labels=FALSE, tck=.02)        ## interior ticks
    }
    

    ## pairs
    if (print.mca) {
      if (is.null(mca)) warning("mca component not found.")
      else {
        yy <- mca$height
        lower <- mca$table[,"lower"]
        upper <- mca$table[,"upper"]
        lower <- ifelse(lower == ( - Inf), par()$usr[1], lower)
        upper <- ifelse(upper == Inf, par()$usr[2], upper)

        signif <- (mca$table[,"lower"] * mca$table[,"upper"] > 0)

        if (any(!signif))
        segments(lower[!signif], yy[!signif]/2,
                 upper[!signif], yy[!signif]/2,
                 xpd=TRUE,
                 lty=lty.mca.not.signif,
                 col=col.mca.not.signif,
                 lwd=lwd.mca.not.signif)

        if (any(signif))
        segments(lower[signif], yy[signif]/2,
                 upper[signif], yy[signif]/2,
                 xpd=TRUE,
                 lty=lty.mca.signif,
                 col=col.mca.signif,
                 lwd=lwd.mca.signif)

        y.01 <- diff(par()$usr[3:4])/100

        if (any(!signif))
        segments(mca$table[!signif,"estimate"], yy[!signif]/2-y.01,
                 mca$table[!signif,"estimate"], yy[!signif]/2+y.01,
                 xpd=TRUE, lty=1, col=col.mca.not.signif)

        if (any(signif))
        segments(mca$table[signif,"estimate"], yy[signif]/2-y.01,
                 mca$table[signif,"estimate"], yy[signif]/2+y.01,
                 xpd=TRUE, lty=1, col=col.mca.signif)

        if (any(!signif))
          mtext(dimnames(mca$table)[[1]][!signif], at=yy[!signif]/2,
                side=4, adj=0,
                line=.5,
                col=col.mca.not.signif, las=1,
                xpd=if.R(s=TRUE, r=NA))

        if (any(signif))
          mtext(dimnames(mca$table)[[1]][signif], at=yy[signif]/2,
                side=4, adj=0,
                line=.5,
                col=col.mca.signif, las=1,
                xpd=if.R(s=TRUE, r=NA))

        axis(4, at=yy/2, tck=-.01, labels=FALSE)
      }
    }

    ## lmat argument
    if (print.lmat) {
      if (is.null(lmat)) warning("lmat component not found.")
      else {
        yy <- lmat$height
        lower <- lmat$table[,"lower"]
        upper <- lmat$table[,"upper"]
        lower <- ifelse(lower == ( - Inf), par()$usr[1], lower)
        upper <- ifelse(upper == Inf, par()$usr[2], upper)
        
        signif <- (lmat$table[,"lower"] * lmat$table[,"upper"] > 0)
        
        if (any(!signif))
        segments(lower[!signif], yy[!signif]/2,
                 upper[!signif], yy[!signif]/2,
                 xpd=TRUE,
                 lty=lty.lmat.not.signif,
                 col=col.lmat.not.signif,
                 lwd=lwd.lmat.not.signif)

        if (any(signif))
        segments(lower[signif], yy[signif]/2,
                 upper[signif], yy[signif]/2,
                 xpd=TRUE,
                 lty=lty.lmat.signif,
                 col=col.lmat.signif,
                 lwd=lwd.lmat.signif)
        
        y.01 <- diff(par()$usr[3:4])/100

        if (any(!signif))
        segments(lmat$table[!signif,"estimate"], yy[!signif]/2-y.01,
                 lmat$table[!signif,"estimate"], yy[!signif]/2+y.01,
                 xpd=TRUE, lty=1, col=col.lmat.not.signif)

        if (any(signif))
        segments(lmat$table[signif,"estimate"], yy[signif]/2-y.01,
                 lmat$table[signif,"estimate"], yy[signif]/2+y.01,
                 xpd=TRUE, lty=1, col=col.lmat.signif)
        
        if (any(!signif))
          mtext(dimnames(lmat$table)[[1]][!signif], at=yy[!signif]/2,
                side=4, adj=1, line=-.6, col=col.lmat.not.signif, las=1)

        if (any(signif))
          mtext(dimnames(lmat$table)[[1]][signif], at=yy[signif]/2,
                side=4, adj=1, line=-.6, col=col.lmat.signif, las=1)
        
        axis(4, at=yy/2, tck=.01, labels=FALSE)
      }
    }
    ## par(old.par) ## restore par() settings
    invisible(c(main.method.phrase, main2.method.phrase))
  }
