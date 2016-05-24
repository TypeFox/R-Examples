plot.gv <-
function(x, line.res=100, pch=1,
         legend=TRUE, leg.x=NA, leg.y=NA, leg.cex=1, ...) {
### TODO: check if it is a multi should be via a class attribut
### or in the object parameters slot.
### plot.gv can hadle multi and single gv objects

    ## check if gv is multi
    multi <- FALSE
    if ("gamma.mat" %in% names(x)) multi <- TRUE

    ## Get all distance classes to X-axis
    X <- x$lag[!is.na(x$lag)] 
    Y <- x$gamma[!is.na(x$gamma)]

    # Check if a model is present
    mtest <- mtest.gv(x)
    if (mtest) {
        xx <- seq(0, ceiling(max(X)), length.out = line.res)
        yy <- predict(x, xx)
    }

    # Prepare point size relative to n
    lab.n <- range(x$n)
    cex <- 2 * x$n / lab.n[2]
    
    ## ### Begin plotting ### ##
    plot.new()

    ## Get X- and Y- axes ranges
    x.range <- c(0, max(X))
    mY <- max(Y)
    if (multi) mY <- max(x$gamma.ci[2,])
    if (mtest) {
        y.range <- c(0, max(c(mY,yy)))
    } else {
        y.range <- c(0, mY)
    }

    plot.window(x.range, y.range)
    if (multi) {
        lines(X, x$gamma.ci[1,], lty=2)
        lines(X, x$gamma.ci[2,], lty=2)
        lines(X, Y, lty=1)
    }
    
    points(X,Y, pch=pch, cex=cex*leg.cex, ...)

    if (mtest) lines(xx, yy, col='red', ...)

    axis(1)
    axis(2)
    box()

    # Plot Legend
    if (legend) {
        leg <- pretty(x$n)
        leg[1] <- 1
        if (is.na(leg.x)) leg.x = 0.9 * diff(x.range)
        if (is.na(leg.y)) leg.y = 0.06 * length(leg) * diff(y.range)
        legend(leg.x, leg.y, legend=leg, pch=pch, 
               pt.cex = leg.cex * 2 * leg / lab.n[2],
               title = expression(italic('n')*' size'))
    }

    # Plot titles
    main <- "Semi-Variogram"
    if (mtest) 
        main <- paste(main, "\n",
                      "Model:", x$model$type,
                      "Sill:", round(x$model$sill, 3),
                      "Range:", round(x$model$range, 3),
                      "Nugget:", round(x$model$nugget, 3))
    title(main=main, xlab='Distance', ylab='Semivariance')
}


