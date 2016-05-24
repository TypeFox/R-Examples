plot.glc <- 
    function(x, fitdb = TRUE, initdb = FALSE, xlim = NULL, ylim = NULL, bg, pch, ...)
{
    if (!inherits(x, "glc")) stop("object not of class \"glc\"")
    data <- x$model
    Terms <- x$terms
    X <- model.matrix(delete.response(Terms), data)
    res <- model.response(data)
    categ <- x$category
    lev <- levels(factor(categ))
    tlabels <- attr(Terms,"term.labels")
    initpar <- x$initpar
    params <- x$par
    names(initpar$coeffs) <- tlabels
    names(params$coeffs) <- tlabels
    xint <- match("(Intercept)", colnames(X), nomatch=0L)
    if(xint > 0) X <- X[, -xint, drop=FALSE]
    if(missing(bg)) bg=c("white","gray")[res]
    if(missing(pch)) pch=c(21,24)[categ]
    
    if(ncol(X) == 1)
    {
        if(is.factor(categ)) categ <- c(1:(length(lev)))[categ]
        if(length("xlim") != 2) xlim=c(.5, 2.5)
        plot(x=categ, y=X, type='p', xaxt='n', xlim=xlim, bg=bg, pch=21, ...)
	axis(1, at=c(1:2), lev)
	if(initdb) abline(h = coef(initpar), col = "blue")
        if(fitdb) abline(h = coef(params), col = "red")
    } else if(ncol(X) == 2) {
        plot(X[,1L], X[,2L], type='p', bg=bg, pch=pch, 
            xlab=tlabels[1L], ylab=tlabels[2L], xlim=xlim, ylim=ylim,...)
        if(initdb) abline(coef(initpar), col="blue")
        if(fitdb) abline(coef(params), col="red")
    } else if(ncol(X) == 3) {
        plot3d.glc(x, fitdb=fitdb, initdb=initdb, ...)
    } else {
	stop("unsupported number of dimentions.")
    }
    invisible(NULL)
}

plot3d.glc <- 
    function(x, 
        fitdb = TRUE, 
        initdb = FALSE, 
        lims = NULL, 
        alpha = .5, ...)
{
    if (!inherits(x, "glc")) stop("object not of class \"glc\"")
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
    #a function to obtain the limit of xyz coordinate
    fz <- function(params, xlim, ylim){
        fzc <- c(params$coeffs[1:2],params$bias)/-params$coeffs[3]
        x <- c(xlim[1],xlim[2],xlim[2],xlim[1])
        y <- c(ylim[1],ylim[1],ylim[2],ylim[2])
        z <- fzc[1] * x + fzc[2] * y + fzc[3]
        cbind(x,y,z)
    }
    if(is.null(lims)) {
        lims <- sapply(as.data.frame(X), function(x){extendrange(range(x, finite=TRUE))})
    } else if (is.list(lims)) {
        lims <- cbind(lims[[1]], lims[[2]], lims[[3]])
    }
    
    pla0 <- fz(initpar, xlim=lims[,1],ylim=lims[,2])
    pla1 <- fz(params, xlim=lims[,1],ylim=lims[,2])
    lims[,3] <- c(min(c(pla0[,3], pla1[,3])), max(pla0[,3], pla1[,3]))
    rgl::open3d()
    rgl::plot3d(X, type="n", xlim=lims[,1], ylim=lims[,2], zlim=lims[,3])
    rgl::points3d(X, col=c("black","gray")[categ])
    if(initdb) rgl::quads3d(pla0, col="red", alpha=alpha)
    if(fitdb) rgl::quads3d(pla1, col="blue", alpha=alpha)
    invisible(NULL)
}


