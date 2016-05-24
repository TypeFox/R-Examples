grid.image <- 
function(intpl, grid, breaks=10, ic=1, colFUN=heat.colors, 
         main=colnames(intpl)[ic], xlab=colnames(grid)[1], 
         ylab=colnames(grid)[2], sclab=NA, ...) {

    x <- sort(unique(grid[,1]))
    y <- sort(unique(grid[,2]))
    res <- c(diff(x[1:2]), diff(y[1:2]))
    rng.x <- range(x)
    rng.y <- range(y)
    
    if (is.vector(intpl)) {
        # vector found, ignore ic value.
        mat <- tapply(intpl, list(grid[,1], grid[,2]), mean)
    } else {
        if (is.numeric(ic) & ic > ncol(intpl)) 
            stop("ic cannot be higher than columns available.")
        if (is.character(ic) & !(ic %in% colnames(intpl)))
            stop("ic does not correspond to any column name.")     
        mat <- tapply(intpl[,ic], list(grid[,1], grid[,2]), mean)
    }

    # plot scale first
    par(fig = c(0,1,0,0.3), ...)
    plot.new()
    plot.window(range(mat, na.rm=TRUE), c(0,1))

    #Automatic scale number of breaks if binary data is found
    test <- na.exclude(unique(as.vector(mat)))
    if (length(test) == 2) breaks <- 1

    sx <- seq(min(mat, na.rm=TRUE), 
              max(mat, na.rm=TRUE), 
              length.out=breaks+1)
    c.len <- diff(range(sx))/breaks

    if (breaks < 1) {
        stop("Breaks must be 1 or higher")
    } else if (breaks == 1) {
        lbl.s <- sx
        corr.pos <- c(-c.len/4, c.len/4)
    } else if (breaks > 5) {
        lbl.s <- sx[round(seq(1, length(sx), length.out=6), 0)]
        corr.pos <- c.len*0.5
    } else {
        lbl.s <- sx
        corr.pos <- c.len*0.5
    }

    image(sx, 1, matrix(seq(0,1,length.out=breaks+1)), 
          col=colFUN(breaks+1), add=T)
    axis(1, at = lbl.s-corr.pos, labels = round(lbl.s, 2), 
         cex.axis = par('cex')*0.75)

    usr <- par('usr')
    if (usr[1] < min(sx) - c.len*0.5) usr[1] = min(sx) - c.len*0.5
    if (usr[2] > max(sx) + c.len*0.5) usr[2] = max(sx) + c.len*0.5
    rect(usr[1], usr[3], usr[2], usr[4])
    axis(1, at = mean(sx), labels = sclab, pos = -1, tick = FALSE, 
         cex.axis = par('cex')*0.75)

    # plot interpolation results
    par(fig = c(0,1,0.2,1), new=TRUE)

    plot.new()
    plot.window(rng.x, rng.y, asp=1)

    # plot matrix
    image(x, y, mat, col=colFUN(breaks+1), add=TRUE) 

    # plot axis
    lbl.x <- pretty(x)
    lbl.y <- pretty(y)
    #remove extra points from labels for coeherent axis
    if (lbl.x[1] < rng.x[1]) lbl.x <- lbl.x[-1]
    if (lbl.x[length(lbl.x)] > rng.x[2]) lbl.x <- lbl.x[-length(lbl.x)]
    if (lbl.y[1] < rng.y[1]) lbl.y <- lbl.y[-1]
    if (lbl.y[length(lbl.y)] > rng.y[2]) lbl.y <- lbl.y[-length(lbl.y)]

    area <- c(rng.x[1]-res[1]/2, rng.y[1]-res[2]/2, rng.x[2]+res[1]/2, rng.y[2]+res[2]/2)
    axis(1, at=lbl.x, labels=lbl.x, pos=area[2])
    axis(2, at=lbl.y, labels=lbl.y, pos=area[1])
    axis(3, at=lbl.x, labels=FALSE, pos=area[4])
    axis(4, at=lbl.y, labels=FALSE, pos=area[3])
    rect(area[1], area[2], area[3], area[4])

    # plot labels
    title(main=main)
    axis(1, at=mean(x), labels=xlab, pos=area[2]-diff(rng.y)*0.1, tick = FALSE)
    axis(2, at=mean(y), labels=ylab, pos=area[1]-diff(rng.x)*0.1, tick = FALSE)


}

