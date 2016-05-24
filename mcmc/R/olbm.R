
olbm <- function(x, batch.length, demean = TRUE) {
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    storage.mode(x) <- "double"
    if (batch.length > n) stop("batch.length must be <= nrow(x)")
    if (demean) {
    	mean <- apply(x, 2, mean)
    	no.calc.mean <- TRUE
    } else {
    	mean <- double(p)
    	no.calc.mean <- FALSE
    }
    out <- .C("olbm",
    	x=x,
    	n=as.integer(n),
    	p=as.integer(p),
    	batch.length=as.integer(batch.length),
    	mean=as.double(mean),
    	var=matrix(as.double(0), p, p),
    	no.calc.mean=as.logical(no.calc.mean), PACKAGE = "mcmc")
    return(out$var)
}

