#generic function
bride <- function(p, y, n.bins=10, ...) UseMethod("bride")

#default S3 method
bride.default <- function(p, y, n.bins=10, ...)
{
    p <- as.vector(p)
    y <- as.vector(y)

    #error checking
    if (length(y) != length(p)) {
	stop("y and p must be of equal length.")
    }
    if ( all(y==0) || all(y==1) || any(!(y %in% c(0,1))) ) {
	stop(paste("y can be only 1 ('event') or 0 ('non-event').",
              "Also, y must contain at least one 1 and one 0."))
    }
    if (n.bins < 2 || n.bins >= length(p)) {
        stop("n.bins must be between 2 and length(p)-1.")
    }

    #calculate
    brd <- CalculateBrierDecomp(p=p, y=y, n.bins=n.bins)
    brd$call <- match.call()
    class(brd) <- "bride"
    brd
}

#print method
print.bride <- function(x, ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\n")
    cat(sprintf("REL = %1.4f +/- %1.4f\n", x$rel2, sqrt(x$rel2.var)))
    cat(sprintf("RES = %1.4f +/- %1.4f\n", x$res2, sqrt(x$res2.var)))
    cat(sprintf("UNC = %1.4f +/- %1.4f\n", x$unc2, sqrt(x$unc2.var)))
    cat("\n")
    cat(sprintf("REL - RES + UNC = %1.4f\n", x$rel2-x$res2+x$unc2))
    cat(sprintf("Br = %1.4f\n", mean( (x$y-x$p)^2 )))
    cat("\n")
}

#summary method
summary.bride <- function(object, ...) {
    print.bride(object)
}

