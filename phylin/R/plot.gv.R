plot.gv <-
function(x, line.res=100, pch=1,
         legend=TRUE, leg.x=NA, leg.y=NA, leg.cex=1, ...) {
    X <- x$lag[!is.na(x$lag)]
    Y <- x$gamma[!is.na(x$gamma)]

    mtest <- mtest.gv(x)
    if (mtest) {
        xx <- seq(0, ceiling(max(X)), length.out = line.res)
        yy <- predict(x, xx)
    }

    #Prepare point size relative to n
    lab.n <- range(x$n)
    cex <- 2 * x$n / lab.n[2]
    
    # Begin plot
    plot.new()
    x.range <- c(0, max(X))
    if (mtest) {
        y.range <- c(0, max(c(Y,yy)))
    } else {
        y.range <- c(0, max(Y))
    }

    plot.window(x.range, y.range)
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


