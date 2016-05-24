image.dist <-
function(x, grad, lab=TRUE, ...) {
    ## labels not yet used but could be done to add labels to left axis
    z <- as.matrix(x)
    make.lab <- !is.logical(lab) || lab
    if (!is.logical(lab)) {
        if (length(lab) != nrow(z))
            stop("inadequate length for 'lab'")
        make.lab <- TRUE
    }
    if (is.logical(lab) && lab) {
        lab <- labels(x)
        make.lab <- TRUE
    }
    if (!missing(grad)) {
        z <- z[order(grad), order(grad)]
        if (make.lab)
            lab <- lab[order(grad)]
    }
    z[upper.tri(z)] <- diag(z)[1]
    pv <- t(1-z/max(z))
    image(pv, ylim=c(1,0), axes=FALSE, ...)
    if (make.lab) {
        op <- par(las=1)
        mtext(lab, 2, at=(1:length(lab)-0.5)/length(lab))
        par(op)
     }
    invisible(pv)
}

