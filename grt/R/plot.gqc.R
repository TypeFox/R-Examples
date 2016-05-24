plot.gqc <- 
    function(x, fitdb=TRUE, initdb=FALSE, xlim=NULL, ylim=NULL, bg, pch, npoints = 100, ...)
{
    if (!inherits(x, "gqc")) stop("object not of class \"gqc\"")
    data <- x$model
    Terms <- x$terms
    X <- model.matrix(delete.response(Terms), data)
    xint <- match("(Intercept)", colnames(X), nomatch=0L)
    if(xint > 0) X <- X[, -xint, drop=FALSE]
    res <- model.response(data)
    categ <- x$category
    lev <- levels(factor(categ))
    labels <- attr(Terms,"term.labels")
    
    if(missing(bg)) bg <- c("white","gray")[res]
    if(missing(pch)) pch <- c(21,24)[categ]
    if(missing(xlim) | is.null(xlim)) 
        xlim <- extendrange(range(X[,1]), f = 0.1)
    if(missing(ylim) | is.null(ylim)) 
        ylim <- extendrange(range(X[,2]), f = 0.1)
    
    if(ncol(X) == 2) {
        plot(X[,1L], X[,2L], type="p", bg=bg, pch=pch, 
            xlab=labels[1L], ylab=labels[2L], xlim=xlim, ylim=ylim,...)
        if(initdb)
            lines(x$initpar, xlim = xlim, ylim = ylim, col = "blue")
        if(fitdb)
            lines(x$par, xlim = xlim, ylim = ylim, col = "red")
    }
    invisible(NULL)
}

lines.gqcStruct <- function(x, 
    xlim = c(0,1), ylim = c(0,1), 
    npoints = 100, col = "black", 
    ...)
{
    if (!inherits(x, "gqcStruct")) stop("object not of class \"gqcStruct\"")
    if (length(x$coeffs) != 5) stop("too many/few coeffs")
    xlim <- extendrange(xlim, f = 0.1)
    ylim <- extendrange(ylim, f = 0.1)
    xseq <- seq(xlim[1], xlim[2], length.out=npoints)
    yseq <- seq(ylim[1], ylim[2], length.out=npoints)
    x1 <- matrix(rep(xseq, npoints), ncol=npoints, byrow=FALSE)
    x2 <- matrix(rep(yseq, npoints), ncol=npoints, byrow=TRUE)
    a <- c(x$coeffs, x$bias)
    z <- a[1]*x1^2 + a[2]*x2^2 + a[3]*x1*x2 + a[4]*x1 + a[5]*x2 + a[6]
    #contour(x = xseq, y = yseq, z = z, levels=0, drawlabels=F, add = T, col = col)
    xy <- contourLines(x = xseq, y = yseq, z = z, levels=0)[[1]][-1]
    lines(xy, col = col, ...)
    invisible(xy)
}

plot3d.gqc <- 
    function(x, fitdb = TRUE, initdb = FALSE, 
        lims = NULL, npoints = 100, alpha = .5, 
        fill = TRUE, smooth = FALSE, ...)
{
    if (!inherits(x, "gqc")) stop("object not of class \"gqc\"")
    data <- x$model
    Terms <- x$terms
    X <- model.matrix(delete.response(Terms), data)
    res <- model.response(data)
    categ <- x$category
    lev <- levels(factor(categ))
    labels <- attr(Terms,"term.labels")
    initpar <- x$initpar
    params <- x$par
    names(initpar$coeffs) <- labels
    names(params$coeffs) <- labels
    xint <- match("(Intercept)", colnames(X), nomatch=0L)
    if(xint > 0) X <- X[, -xint, drop=FALSE]
    if(ncol(X) != 3) stop("ncol(X) is not 3")
    if(is.null(lims)) {
        lims <- sapply(as.data.frame(X), function(x){extendrange(range(x, finite=TRUE), f = 0.1)})
    } else if (is.list(lims)) {
        lims <- cbind(lims[[1]], lims[[2]], lims[[3]])
    }
    
    rgl::plot3d(X, type="n", xlim=lims[,1], ylim=lims[,2], zlim=lims[,3])
    rgl::points3d(X, col=c("black","gray")[categ])

    #plotting the 3d surface(s)
    x1 <- seq(lims[1,1], lims[2,1], length.out=npoints)
    x2 <- seq(lims[1,2], lims[2,2], length.out=npoints)
    x3 <- seq(lims[1,3], lims[2,3], length.out=npoints)
    if(initdb){
        a <- c(initpar$coeffs, initpar$bias)
        f <- function(x, y, z){
             a[1]*x^2 + a[2]*y^2 + a[3]*z^2 + a[4]*x*y + a[5]*x*z + a[6]*y*z + a[7]*x + a[8]*y + a[9]*z + a[10]
        }
        misc3d::contour3d(f, level = 0, x = x1, y = x2, z = x3, color = "red", alpha = alpha, 
            fill = fill, material = "default", 
            smooth = smooth, add = TRUE, draw = TRUE, engine = "rgl", 
            xlim=lims[,1], ylim=lims[,2], zlim=lims[,3])
    }
    if(fitdb){
        a <- c(params$coeffs, params$bias)
        f <- function(x, y, z){
             a[1]*x^2 + a[2]*y^2 + a[3]*z^2 + a[4]*x*y + a[5]*x*z + a[6]*y*z + a[7]*x + a[8]*y + a[9]*z + a[10]
        }
        misc3d::contour3d(f, level = 0, x = x1, y = x2, z = x3, color = "blue", alpha = alpha, 
            fill = fill, material = "default", 
            smooth = smooth, add = TRUE, draw = TRUE, engine = "rgl",
            xlim=lims[,1], ylim=lims[,2], zlim=lims[,3])
    }
    invisible(NULL)
}
