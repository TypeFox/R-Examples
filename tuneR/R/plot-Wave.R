plot.Wave.channel <- 
function(x, xunit, ylim, xlab, ylab, main, nr, simplify, axes = TRUE, yaxt = par("yaxt"), 
        las = 1, center = TRUE, ...){
    channel <- if(is(x, "WaveMC")) x@.Data[,1] else x@left
    null <- if(x@bit == 8) 127 else 0
    l <- length(channel)
    if(simplify && (l > nr)){
        nr <- ceiling(l / round(l / nr))
        index <- seq(1, l, length = nr)
        if(xunit == "time") index <- index / x@samp.rate
        mat <- matrix(c(channel, if(l %% nr > 0) rep(NA, nr - (l %% nr))), 
            ncol = nr)
        rg <- apply(mat, 2, range, na.rm = TRUE)
        plot(rep(index, 2), c(rg[1,], rg[2,]), type = "n", yaxt = "n", ylim = ylim, 
            xlab = xlab, ylab = NA, main = main, axes = axes, las = las, ...)
        segments(index, rg[1,], index, rg[2,], ...)
    }
    else{
        index <- seq(along = channel)
        if(xunit == "time") index <- index / x@samp.rate
        plot(index, channel,
            type = "l", yaxt = "n", ylim = ylim, xlab = xlab, 
            ylab = NA, main = main, axes = axes, las = las, ...)
    }
    mtext(ylab, side = 4, line = 0.5, at = mean(par("usr")[3:4]), cex = par("cex.lab"))

    if(!center || all(ylim <= 0)) {
      at <- axTicks(2)
    } else {
      at <- round((ylim[2] - null) * 2/3, -floor(log(ylim[2], 10)))
      at <- null + c(-at, 0, at)
    }
    if(axes) axis(2, at = at, yaxt = yaxt, las = las)
}


setMethod("plot", signature(x = "Wave", y = "missing"),
function(x, info = FALSE, xunit = c("time", "samples"), 
    ylim = NULL, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, 
    simplify = TRUE, nr = 2500, axes = TRUE, yaxt = par("yaxt"), las = 1, 
    center = TRUE, ...){
    
    xunit <- match.arg(xunit)
    if(is.null(xlab)) xlab <- xunit
    stereo <- x@stereo
    l <- length(x@left)
    if(center && is.null(ylim)){
        ylim <- range(x@left, x@right)
        if(x@bit == 8)
            ylim <- c(-1, 1) * max(abs(ylim - 127)) + 127
        else
            ylim <- c(-1, 1) * max(abs(ylim))
    }
    if(stereo){
        if(length(ylab)==1) ylab <- rep(ylab, 2)
        opar <- par(mfrow = c(2,1), 
            oma = c((if(info) 6.1 else 5.1) + if(!is.null(sub)) 0.5 else 0, 0, 4.1, 0))
        on.exit(par(opar))
        mar <- par("mar")
        par(mar = c(0, mar[2], 0, mar[4]))
        plot.Wave.channel(mono(x, "left"), xunit = xunit,
            ylab = if(is.null(ylab)) "left channel" else ylab[1], 
            main = NULL, sub = NULL, xlab = NA, ylim = ylim, 
            xaxt = "n", simplify = simplify, nr = nr, 
            axes = axes, yaxt = yaxt, las = las, center = center, ...)
        plot.Wave.channel(mono(x, "right"), xunit = xunit,
            ylab = if(is.null(ylab)) "right channel" else ylab[2],
            main = NULL, sub = sub, xlab = NA, ylim = ylim,  
            simplify = simplify, nr = nr, 
            axes = axes, yaxt = yaxt, las = las, center = center, ...)
        title(main = main, outer = TRUE, line = 2)
        title(xlab = xlab, outer = TRUE, line = 3)
        title(sub  = sub , outer = TRUE, line = 4)
        par(mar = mar)
    }
    else{
        if(info){
            opar <- par(oma = c(2, 0, 0, 0))
            on.exit(par(opar))
        }
        plot.Wave.channel(x, xunit = xunit, 
            ylab = if(is.null(ylab)) "" else ylab,
            main = main, sub = sub, xlab = xlab, ylim = ylim,
            simplify = simplify, nr = nr, 
            axes = axes, yaxt = yaxt, las = las, center = center, ...)            
    }
    if(info){
        mtext(paste("Wave Object: ",  
                l, " samples (", 
                round(l / x@samp.rate, 2),  " sec.), ",
                x@samp.rate, " Hertz, ",
                x@bit, " bit, ",
                if(stereo) "stereo." else "mono.", sep = ""), 
            side = 1, outer = TRUE, line = (if(stereo) 5 else 0) + if(!is.null(sub)) 0.5 else 0, ...)
    }
})


setMethod("plot", signature(x = "WaveMC", y = "missing"),
function(x, info = FALSE, xunit = c("time", "samples"), 
    ylim = NULL, main = NULL, sub = NULL, xlab = NULL, ylab = colnames(x), 
    simplify = TRUE, nr = 2500, axes = TRUE, yaxt = par("yaxt"), las = 1, 
    center = TRUE, mfrow = NULL, ...){
    
    xunit <- match.arg(xunit)
    if(is.null(xlab)) xlab <- xunit
    l <- nrow(x)
    if(center && is.null(ylim)){
        ylim <- range(x)
        if(x@bit == 8)
            ylim <- c(-1, 1) * max(abs(ylim - 127)) + 127
        else
            ylim <- c(-1, 1) * max(abs(ylim))
    }
    if(is.null(mfrow)) {
      mfrow <- switch(as.character(ncol(x)),
        "2" = c(2,1),
        "3" = c(3,1),
        "4" = c(4,1),
        "5" = c(5,1),
        "6" = c(3,2),
        "7" = c(4,2),
        "8" = c(4,2),
        "9" = c(3,3),
        "10" = c(4,3),
        "11" = c(4,3),
        "12" = c(4,3),
        "13" = c(5,3),
        "14" = c(5,3),
        "15" = c(5,3),
        "16" = c(4,4),
        c(5,4)
      )
    }
      
    if(ncol(x)==1) {
      if(info){
              opar <- par(oma = c(2, 0, 0, 0))
              on.exit(par(opar))
      }
      plot.Wave.channel(x, xunit = xunit, 
          ylab = if(is.null(ylab)) "" else ylab,
          main = main, sub = sub, xlab = xlab, ylim = ylim,
          simplify = simplify, nr = nr, 
          axes = axes, yaxt = yaxt, las = las, center = center, ...)
    } else {
      opar <- par(mfrow = mfrow, oma = c((if(info) 6.1 else 5.1) + if(!is.null(sub)) 0.5 else 0, 0, 4.1, 0))
      on.exit(par(opar))
      if(length(ylab) == 1) ylab <- rep(ylab, ncol(x))
      mar <- par("mar")
      par(mar = c(0, mar[2], 0, mar[4]))
      for(i in 1:ncol(x)) {
        plot.Wave.channel(x[,i], xunit = xunit,
            ylab = if(is.null(ylab)) paste("channel ", i) else ylab[i],         
            main = NULL, sub = NULL, xlab = NA, ylim = ylim,   
            xaxt = if(i %in% (ncol(x) - 0:(mfrow[2]-1))) "s" else "n", simplify = simplify, nr = nr, 
            axes = axes, yaxt = yaxt, las = las, center = center, ...)

        if(i %in% (ncol(x) - 0:(mfrow[2]-1))) {
            oparxpd <- par(xpd=NA)
            title(xlab=xlab)
            par(oparxpd)
        }
      }

      title(main = main, outer = TRUE, line = 2)
      title(sub  = sub , outer = TRUE, line = 4)
      par(mar = mar)
    }
    if(info){
        mtext(paste("WaveMC Object: ",  
                l, " samples (", 
                round(l / x@samp.rate, 2),  " sec.), ",
                x@samp.rate, " Hertz, ",
                x@bit, " bit, ",
                ncol(x), " Channels.", sep=""), 
            side = 1, outer = TRUE, line = (if(ncol(x)==1) 0 else 5) + if(!is.null(sub)) 0.5 else 0, ...)
    }
})
