
# stuff for improving the presentation of tables, etc.
# a.k.a. bertin matrices.
#
# (C) ceeboo 2005, 2006

# the interface to the stress functions allows for
# arbitrary subsetting (see the wrapper in C).

stress <- function(x, rows=NULL, cols=NULL, type="moore") {
    TYPE <- c(1,2,3)
    names(TYPE) <- c("moore", "neumann")
    if (inherits(x, "dist"))
        x <- as.matrix(x)
    if (!is.matrix(x))
       stop(paste(sQuote("x"),"not a matrix"))
    if (!is.double(x))
       storage.mode(x) <- "double"
    if (is.null(rows))
       rows <- as.integer(1:dim(x)[1])
    if (is.null(cols))
       cols <- as.integer(1:dim(x)[2])
    type <- as.integer(TYPE[type])
    x <- .Call(R_stress, x, rows, cols, type)
    x
}

# interface to distance computation based on the above
# stress functions (auto-distances only)

stress.dist <- function(x, rows=NULL, cols=NULL, bycol=FALSE, type="moore") {
    TYPE <- c(1,2)
    names(TYPE) <- c("moore", "neumann")
    if (inherits(x, "dist"))
        as.matrix(x)
    if (!is.matrix(x))
       stop(paste(sQuote("x"),"not a matrix"))
    if (!is.double(x))
       storage.mode(x) <- "double"
    if (is.null(rows))
       rows <- as.integer(1:dim(x)[1])
    if (is.null(cols))
       cols <- as.integer(1:dim(x)[2])
    type <- as.integer(TYPE[type])
    storage.mode(bycol) <- "logical"
    #
    obj <- .Call(R_stress_dist, x, rows, cols, bycol, type)
    # return dist object
    if (bycol)
    obj <- structure(obj, Size= if (bycol) dim(x)[2] else dim(x)[1], 
		             class="dist", Diag=FALSE, Upper=FALSE, 
		             Labels= if (bycol) { if (is.null(colnames(x))) cols 
                                          else  colnames(x)[cols] } 
                             else { if (is.null(rownames(x))) rows 
                                    else rownames(x)[rows] },
		             method=names(TYPE[type]))
    obj
}

# reorder table like objects (we may use S3 dispatch in the 
# future

order.dist <- function(x, index = FALSE) {
    if (!inherits(x, "dist"))
        stop("'x' not of class dist")
    k <- .Call(R_orderTSP, x, sample(attr(x, "Size")))
    cat("length:", order.length(x, k),"\n")
    if (index)
        return(k)
    subset(x, k)
}

order.matrix <-
function(x, type = "neumann", by = c("both","rows","cols"), index = FALSE) {
    if (!is.matrix(x))
        stop("'x' not a matrix")
    by <- match.arg(by)
    if (by == "both") {
        r <- sample(dim(x)[1])
        c <- sample(dim(x)[2])
        c <- c[.Call(R_orderTSP, stress.dist(x,r,c,TRUE, type), seq(c))] 
        r <- r[.Call(R_orderTSP, stress.dist(x,r,c,FALSE,type), seq(r))]
    } else
    if (by == "rows") {
        r <- sample(dim(x)[1])
        c <- seq(dim(x)[2])
        r <- r[.Call(R_orderTSP, stress.dist(x,r,c,FALSE,type), seq(r))]
    } else
    if (by == "cols") {
        r <- seq(dim(x)[1])
        c <- sample(dim(x)[2])
        c <- c[.Call(R_orderTSP, stress.dist(x,r,c,TRUE, type), seq(c))] 
    }
    cat("stress:",stress(x,r,c,type),"\n")
    if (index)
        return(list(rows=r, cols=c))
    x <- x[r,c]
    if (is.null(rownames(x)))
        rownames(x) <- r
    if (is.null(colnames(x)))
        colnames(x) <- c
    x
}

order.data.frame <-
function(x, type = "neumann", by = c("both","rows","cols"), index = FALSE) {
    if (!inherits(x, "data.frame"))
        stop("'x' not a data frame")
    k <- sapply(x, function(x) is.numeric(x) || is.logical(x))
    if (!any(k)) {
        warning("cannot order on ordinal attributes only")
        if (index)
            return(list(rows=seq(dim(x)[1]),cols=seq(dim(x)[2])))
        x
    }
    z <- as.matrix(as.data.frame(lapply(x[k], function(x) {
        if (is.logical(x))
            as.integer(x)
        else {
            m <- min(x)
            (x+m)/(max(x)-m)
        }
    })))
    o <- order.matrix(z, type, by, index=TRUE)
    if (by == "cols" || by == "both") {
        c <- o$cols
        o$cols <- seq(k)
        o$cols[k] <- c
    }
    if (index)
        return(o)
    x[o$rows,o$cols]
}

### the end
