
# wrapper functions for the Rock algorithm.
#
# note that the behavior for other than the binary distance functions
# has not been tested. therefore, the default relationship between beta 
# and theta may not be meaningful in all cases. 
#
# (C) ceeboo 2005

# compute link counts (internal function)
#
# let me stress that the semantics are unscaled 
# similarities but we package as a dist object 
# for possible future use in different contexts.

rockLink <- function(x, beta=0.5) {
    if (!inherits(x, "dist"))
       stop(paste(sQuote("x"),"not of class dist"))
    if (!is.double(x))
       storage.mode(x) <- "double"
    storage.mode(beta) <- "double"
    obj <- .Call(R_rockLink, x, beta)
    obj <- structure(obj, Size=attr(x,"Size"),
                     class="dist", Diag=FALSE, Upper=FALSE,
                     Labels=attr(x, "Labels"), method="rock")
    #invisible(obj)
    obj
}

# merge into clusters (internal function)

rockMerge <- function(x, n, theta=0.5, debug=FALSE) {
    if (!inherits(x, "dist"))
       stop(paste(sQuote("x"),"not of class dist"))
    if (n < 1)
       stop(paste(sQuote("n"),"illegal value"))
    if (theta < 0 || theta >= 1)
       stop(paste(sQuote("theta"),"illegal value"))
    if (!is.integer(x))
       storage.mode(x) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(theta) <- "double"
    storage.mode(debug) <- "logical"
    obj <- .Call(R_rockMerge, x, n, theta, debug)
    names(obj) <- c("cl","size")
    names(obj$cl) <- attr(x,"Labels")
    invisible(obj)
}

# classify based on distances to clustered samples 
# (we have to compute these separately; for an
# example wrapper see below; internal function)

rockClass <- function(x, cl, beta=1-theta, theta=0.5) {
    if (!is.matrix(x))
       stop(paste(sQuote("x"),"not a mtrix"))
    if (!is.factor(cl))
       stop(paste(sQuote("cl"),"not a factor"))
    if (!is.double(x))
       storage.mode(x) <- "double"
    storage.mode(beta) <- storage.mode(theta) <- "double"
    storage.mode(cl) <- "integer"
    obj <- .Call(R_rockClass, x, cl, beta, theta)
    names(obj) <- c("cl","size")
    names(obj$cl) <- rownames(x)
    invisible(obj)
}

# cluster interface

rockCluster <- function(x, n, beta=1-theta, theta=0.5, fun="dist", 
                        funArgs=list(method="binary"), debug=FALSE) {
    if (!is.matrix(x))
       warning(paste(sQuote("x"),"not a matrix"))
    if (n < 1)
       stop(paste(sQuote("n"),"illegal value"))
    if (is.function(fun))
       fun <- deparse(substitute(fun))
    # cluster
    cat("Clustering:\n")
    cat("computing distances ...\n")
    rc <- do.call(fun, c(list(x=x), as.list(funArgs)))
    cat("computing links ...\n")
    rc <- rockLink(rc, beta)
    cat("computing clusters ...\n")
    rc <- rockMerge(rc, n, theta, debug)
    rc <- list(x=x, cl=rc$cl, size=rc$size,
               beta=beta, theta=theta, fun=fun, funArgs=funArgs)
    class(rc) <- "rock"
    rc
}

# wrapper for predicting the class of new (or existing) samples
#

predict.rock <- function(object, x, drop=1, ...) {
    if (!is.matrix(x))
       warning(paste(sQuote("x"),"not a matrix"))
    # drop 
    if (drop > 0) {
       d <- which(object$size <= drop)
       if (length(d) > 0) {
          cat("dropping",length(d),"clusters\n")
          object$size <- object$size[-d]
          k <- !object$cl %in% d            # keep
          object$cl <- factor(object$cl[k]) # enforce contiguous indexing !!!
          object$x <- object$x[k,]
       }
    }
    # classify
    cat("computing distances ...\n")
    x <- do.call(object$fun, c(list(x=x, y=object$x), as.list(object$funArgs)))
    cat("computing classes ...\n")
    x <- rockClass(x, object$cl, object$beta, object$theta)
    x
}

fitted.rock <- function(object, ...)
    predict.rock(object, object$x)

print.rock <- function(x, ...) {
    cat(" data:",dim(x$y)[1],"x",dim(x$y)[2],"\n")
    cat(" beta:",x$beta,"\n")
    cat("theta:",x$theta,"\n")
    cat("  fun:",x$fun,"\n")
    cat(" args:",deparse(x$funArgs, control=NULL),"\n")
    print(x$size)
    invisible(x)
}

### the end

