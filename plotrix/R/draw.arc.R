draw.arc <- function(x=1, y=NULL, radius=1, angle1=deg1*pi/180, angle2=deg2*pi/180,
    deg1=0, deg2=45, n=0.05, col=NA, lwd=NA, ...) 
    {
    if (all(is.na(col)))
        col <- par("col")
    if (all(is.na(lwd)))
        lwd <- par("lwd")
    xylim<-par("usr")
    ymult <- getYmult()
    devunits <- dev.size("px")
    draw.arc.0 <- function(x, y, radius, angle1, angle2, n, col, lwd, ...)
        {
        delta.angle <- (angle2 - angle1)
        if (n != as.integer(n))
            n <- as.integer(1+delta.angle/n) # Divide total angle by desired segment angle to get number of segments
        delta.angle <- delta.angle/n
        angleS <- angle1 + seq(0, length=n) * delta.angle
        angleE <- c(angleS[-1], angle2)
        # Move segment starts/ends so that segments overlap enough to make wide segments
        # not have an open slice in them.  The slice is open by delta.angle*half.lwd.user.
        # That subtends an angle of that/(radius+half.lwd.user) radians, from center.
        # Move segment endpoints by half of that, so together they equal that.
        if (n > 1)
            {
            half.lwd.user <- (lwd/2)*(xylim[2]-xylim[1])/devunits[1]
            adj.angle = delta.angle*half.lwd.user/(2*(radius+half.lwd.user))
            angleS[2:n] = angleS[2:n] - adj.angle
            angleE[1:(n-1)] = angleE[1:(n-1)] + adj.angle
            }
        p1x <- x + radius * cos(angleS)
        p1y <- y + radius * sin(angleS) * ymult
        p2x <- x + radius * cos(angleE)
        p2y <- y + radius * sin(angleE) * ymult
        segments(p1x, p1y, p2x, p2y, col=col, lwd=lwd, ...)
        }
    xy <- xy.coords(x, y)
    x <- xy$x
    y <- xy$y
    a1 <- pmin(angle1, angle2)
    a2 <- pmax(angle1, angle2)
    angle1 <- a1
    angle2 <- a2
    args <- data.frame(x, y, radius, angle1, angle2, n, col, lwd, stringsAsFactors=FALSE)
    for (i in 1:nrow(args))
        do.call("draw.arc.0", c(args[i, ], ...))
    invisible(args)
    }
