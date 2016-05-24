plotCDF2 <- function(mle, bound, col = gray(seq(0.9, 0.3, len = 30)), 
                     surv=FALSE, key=TRUE, n.key=10, round.key=2, cex.key=0.6, 
                     xlim=NULL, ylim=NULL, zlim=NULL, breaks=NULL, xlab="", 
                     ylab="", main="", sub=""){
    p <- mle$p
    rects <- mle$rects
    m <- length(p)
    ncol <- length(col)

    # Check input
    if (!is.vector(p)) stop("invalid argument 'p'")
    if (sum(p) > 1+1e-5 | sum(p) < 1-1e-5) warning("sum(p) is not equal to 1")
    if (!is.matrix(rects) | ncol(rects) != 4) {
        stop("invalid argument 'rects': it must be a matrix with 4 columns")
    }
    if (sum(rects[, 1] > rects[, 2]) + sum(rects[, 3] > rects[, 4]) >= 1) {
        stop("invalid argument 'rects': each row x1,x2,y1,y2 must satisfy 
x1<=x2 and y1<=y2")
    }
    if (length(p) != nrow(rects)) stop("length(p) must equal nrow(rects)")
    if (!is.null(breaks)) {
        if (length(breaks) != ncol + 1) {
            stop("length(breaks) must equal length(col)+1")
        }
    }
    if (!is.logical(key)) stop("argument 'key' must be logical")
    if (ncol<2) stop("length(col) must be >=2")

    if (is.null(zlim)) zlim <- c(0,sum(p))

    if (bound == "u"){
        x <- rects[,1]
        y <- rects[,3]
    }
    if (bound == "l"){
        x <- rects[,2]
        y <- rects[,4]
    }
    if (is.null(xlim)) {
        xrange <- range(x)
        xlim <- c(xrange[1]-.04*(xrange[2]-xrange[1]), 
                  xrange[2]+.04*(xrange[2]-xrange[1]))
    } else if (xlim[1] >= min(x) | xlim[2] <= max(x)) {
        stop("xlim must properly contain the range of x-coordinates of the mass points")
    }
    if (is.null(ylim)) {
        yrange <- range(y)
        ylim <- c(yrange[1]-.04*(yrange[2]-yrange[1]),
                  yrange[2]+.04*(yrange[2]-yrange[1]))
    }else if (ylim[1] >= min(y) | ylim[2] <= max(y)) {
        stop("ylim must properly contain the range of y-coordinates of the mass points")
    }

    xord <- order(x)
    x <- x[xord] 
    y <- y[xord]
    p <- p[xord] 
    yord <- order(y)
            
    CDF <- matrix(0,nrow=m+1, ncol=m+1)
    for (i in 1:m){
        CDF[,i+1] <- CDF[,i] + p[i]*c(rep(0,yord[i]),rep(1,m-yord[i]+1))
    }

    # Remove duplicate x and y values
    x <- sort(c(xlim, x))
    y <- sort(c(ylim, y))
    dimx <- length(x)-1
    dimy <- length(y)-1
    remove.x <- rep(0, dimx)
    remove.y <- rep(0, dimy)
    for (i in 1:(dimx - 1)) {
        if (x[i + 1] == x[i]) remove.x[i] <- 1
    }
    for (i in 1:(dimy - 1)){
        if (y[i + 1] == y[i]) remove.y[i] <- 1
    }
    x <- x[!remove.x]
    y <- y[!remove.y]
    CDF <- t(CDF)[!remove.x, !remove.y]
    dimx <- nrow(CDF)
    dimy <- ncol(CDF)

    # Create image plot
    if (is.null(breaks)){
        CDFmin <- min(CDF)
        CDFmax <- max(CDF)
        if (!is.null(zlim)){
            CDFmin <- max(CDFmin,zlim[1])
            CDFmax <- min(CDFmax,zlim[2])
        }
        a <- (CDFmax - CDFmin)/(ncol-1)
        breaks <- seq(from=CDFmin-a/2, by=a, to=CDFmax+a/2)
    } 
    if (surv){
       image(x, y, 1-CDF, col = col, zlim=zlim, breaks=breaks,
             xlab = xlab, ylab = ylab, main = main, sub = sub)
    }else{
       image(x, y, CDF, col = col, zlim=zlim, breaks=breaks,
             xlab = xlab, ylab = ylab, main = main, sub = sub)
    }
    box()           

    # Add color key
    if (key) {
        opar <- par(xpd = TRUE)
        xmin = x[1]
        xmax = x[length(x)]
        ymin = y[1]
        ymax = y[length(y)]
        key.x1 <- xmax + 0.01 * (xmax - xmin)
        key.x2 <- xmax + 0.04 * (xmax - xmin)
        key.y <- ymin + c(0:ncol) * (ymax - ymin)/ncol
        label.x <- xmax + 0.025 * (xmax - xmin)
        rect(rep(key.x1, ncol), key.y[1:ncol], rep(key.x2, ncol), 
            key.y[2:(ncol + 1)], col = col, border = NA)
        rect(key.x1, ymin, key.x2, ymax, border = "black")
        if (ncol <= 1.5 * n.key) {
            ind.label <- c(1:(ncol))
        } else {
            a <- ceiling(log2(ncol/(1.5 * n.key)))
            ind.label <- seq(from = 2^a, by = 2^a, to = ncol+1)
        }
        text(label.x, key.y[ind.label], 
             round(breaks, round.key)[ind.label], cex = cex.key)
        par(opar)
    }
}
