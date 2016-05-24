plot.gcjc <- 
    function(x, fitdb = TRUE, initdb = FALSE, xlim = NULL, ylim = NULL, bg, pch, ...)
{
    if (!inherits(x, "gcjc")) stop("object not of class \"gcjc\"")
    data <- x$model
    Terms <- x$terms
    X <- model.matrix(delete.response(Terms), data)
    res <- as.factor( model.response(data) )
    categ <- as.factor(x$category)
    lev <- levels(factor(categ))
    tlabels <- attr(Terms,"term.labels")
    initpar <- x$initpar
    params <- x$par
    xint <- match("(Intercept)", colnames(X), nomatch=0L)
    if(xint > 0) X <- X[, -xint, drop=FALSE]
    if(missing(bg)) bg=c("white","gray")[res]
    if(missing(pch)) pch=c(21,24)[categ]
    
    .coef <- function(obj)
    {
        sum(obj$bias * obj$coeffs * -1)
    }

    plot(X[,1L], X[,2L], type='p', bg=bg, pch=pch, 
            xlab=tlabels[1L], ylab=tlabels[2L], xlim=xlim, ylim=ylim,...)
    
    if(initdb)
    {
        abline(v=.coef(initpar[[1]]), col="blue")
        abline(h=.coef(initpar[[2]]), col="blue")
    }
    if(fitdb)
    {
        abline(v=.coef(params[[1]]), col="red")
        abline(h=.coef(params[[2]]), col="red")
    }
    invisible(NULL)
}

