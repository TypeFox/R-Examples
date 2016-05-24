#
# k-means based on conjugate convex functions using sparse data 
# structures and centering (and optionally standardizing).
#
# for details see the C source code.
#
# (C) ceeboo 2005, 2007

ccfkms_sample <- function(x, n) {
    if (inherits(x, "dgCMatrix"))
        as(t(x[,sample(dim(x)[2],n)]), "matrix")
    else
        x[sample(dim(x)[1],n),]
}

ccfkms <- function (x, n, p=NULL, par=2, max.iter=100, opt.std=FALSE, 
					                        opt.retry=0, debug=FALSE) {
    ## dgRMatrix is currently broken
    if (inherits(x, "dgTMatrix"))
        x <- t(as(x, "dgCMatrix"))
    else
    if (inherits(x, "dgCMatrix"))
        x <- t(x)
    else
    if (!is.matrix(x))
       stop(paste(sQuote("x"), "invalid argument"))
    if (!missing(n) && length(n) != 1)
        stop(paste(sQuote("n"), "invalid argument"))
    if (is.null(p))
       p <- ccfkms_sample(x, n)
    else
       if (!is.matrix(p) || ifelse(inherits(x,"dgCMatrix"),
                                   dim(x)[1], dim(x)[2]) != dim(p)[2])
          stop(paste(sQuote("p"), "invalid argument"))
    if (is.matrix(x) && !is.double(x))
       storage.mode(x) <- "double"
    if (!is.double(p))
       storage.mode(p) <- "double"
    storage.mode(par) <- "double"
    storage.mode(max.iter) <- "integer"
    storage.mode(opt.std) <- storage.mode(debug) <- "logical"
  
    obj <- .Call(R_ccfkms, x, p, par, max.iter, opt.std, debug)
    if (opt.retry > 0) {
       for (i in 1:opt.retry) {
           p <- ccfkms_sample(x,n)
           robj <- .Call(R_ccfkms, x, p, par, max.iter, opt.std, debug)
	       if (robj[[4]] < obj[[4]])
	          obj <- robj
       }
    }
    names(obj) <- c("centers", "size", "cl", "inv.inf")
    rownames(obj$centers) <- names(obj$size) <- levels(obj$cl)
    colnames(obj$centers) <- if (inherits(x, "dgCMatrix")) 
                 rownames(x) else colnames(x)
    names(obj$cl) <- if (inherits(x, "dgCMatrix")) 
        colnames(x)  else rownames(x)
    obj <- c(obj, par=par, opt.std=opt.std)
    class(obj) <- "ccfkms"
    obj
}

predict.ccfkms <- function(object, x, drop=1, ...) {
    if (inherits(x, "dgTMatrix"))
        x <- t(as(x, "dgCMatrix"))
    else
    if (inherits(x, "dgCMatrix"))
        x <- t(x)
    else
    if (!is.matrix(x))
       stop(paste(sQuote("x"), "invalid argument"))
    if (ifelse(inherits(x, "dgCMatrix"), dim(x)[1], dim(x)[2])
        != dim(object$centers)[2])
       stop(paste(sQuote("x"), "and", sQuote("object"), "do not conform"))
    if (drop > 0) {
        d <- which(object$size <= drop)
        if (length(d) > 0) {
            cat("dropping", length(d), "clusters\n")
            object$size <- object$size[-d]
            k <- !object$cl %in% d
            object$cl <- factor(object$cl[k])
        }
    }
    x <- ccfkms(x, p=object$centers, par=object$par, opt.std=object$opt.std, 
                                                     max.iter=1)
    x$par <- x$opt.std <- NULL  # prohibit reuse
    x
}

###
