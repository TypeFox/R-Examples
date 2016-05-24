cmeans <-
function(x, centers, iter.max = 100, verbose = FALSE,
         dist = "euclidean", method = "cmeans", m = 2,
         rate.par = NULL, weights = 1, control = list())
{

    x <- as.matrix(x)
    xrows <- nrow(x)
    xcols <- ncol(x)

    if(missing(centers))
        stop("Argument 'centers' must be a number or a matrix.")

    dist <- pmatch(dist, c("euclidean", "manhattan"))
    if(is.na(dist)) 
        stop("invalid distance")
    if(dist == -1) 
        stop("ambiguous distance")

    method <- pmatch(method, c("cmeans", "ufcl"))
    if(is.na(method)) 
        stop("invalid clustering method")
    if(method == -1) 
        stop("ambiguous clustering method")

    if(length(centers) == 1) {
        ncenters <- centers
        centers <- x[sample(1 : xrows, ncenters), , drop = FALSE]
        if(any(duplicated(centers))) {
            cn <- unique(x)
            mm <- nrow(cn)
            if(mm < ncenters) 
                stop("More cluster centers than distinct data points.")
            centers <- cn[sample(1 : mm, ncenters), , drop = FALSE]
        }
    }
    else {
        centers <- as.matrix(centers)
        if(any(duplicated(centers))) 
            stop("Initial centers are not distinct.")
        cn <- NULL
        ncenters <- nrow(centers)
        if (xrows < ncenters)
            stop("More cluster centers than data points.")
    }

    if(xcols != ncol(centers))
        stop("Must have same number of columns in 'x' and 'centers'.")

    if(iter.max < 1) 
        stop("Argument 'iter.max' must be positive.")

    if(method == 2) {
        if(missing(rate.par)) {
            rate.par <- 0.3
        }
    }

    reltol <- control$reltol
    if(is.null(reltol))
        reltol <- sqrt(.Machine$double.eps)
    if(reltol <= 0)
        stop("Control parameter 'reltol' must be positive.")

    if(any(weights < 0))
        stop("Argument 'weights' has negative elements.")
    if(!any(weights > 0))
        stop("Argument 'weights' has no positive elements.")
    weights <- rep(weights, length = xrows)
    weights <- weights / sum(weights)
    
    ## <FIXME>        
    ## Do we really want to do this?   
    perm <- sample(xrows)
    x <- x[perm, ]
    weights <- weights[perm]
    ## </FIXME>

    initcenters <- centers
    pos <- as.factor(1 : ncenters)
    rownames(centers) <- pos

    if(method == 1) {
        retval <- .C("cmeans",
                     as.double(x),
                     as.integer(xrows),
                     as.integer(xcols),
                     centers = as.double(centers),
                     as.integer(ncenters),
                     as.double(weights),
                     as.double(m),
                     as.integer(dist - 1),
                     as.integer(iter.max),
                     as.double(reltol),
                     as.integer(verbose),
                     u = double(xrows * ncenters),
                     ermin = double(1),
                     iter = integer(1),
                     PACKAGE = "e1071")
    }
    else if(method == 2) {
        retval <- .C("ufcl",
                     x = as.double(x),                         
                     as.integer(xrows),
                     as.integer(xcols), 
                     centers = as.double(centers),
                     as.integer(ncenters),
                     as.double(weights),
                     as.double(m),
                     as.integer(dist - 1),
                     as.integer(iter.max),
                     as.double(reltol),                         
                     as.integer(verbose),
                     as.double(rate.par),
                     u = double(xrows * ncenters),
                     ermin = double(1),
                     iter = integer(1),
                     PACKAGE = "e1071")
    }
        
    centers <- matrix(retval$centers, ncol = xcols,
                      dimnames = list(1 : ncenters,
                                      colnames(initcenters)))
    u <- matrix(retval$u, ncol = ncenters,
                dimnames = list(rownames(x), 1 : ncenters))
    u <- u[order(perm), ]
    iter <- retval$iter - 1
    withinerror <- retval$ermin

    cluster <- apply(u, 1, which.max)
    clustersize <- as.integer(table(cluster))
  
    retval <- list(centers = centers, size = clustersize,
                   cluster = cluster, membership = u, iter = iter,
                   withinerror = withinerror, call = match.call())
    class(retval) <- c("fclust")
    return(retval)
}

print.fclust <-
function(x, ...)
{
    cat("Fuzzy c-means clustering with", length(x$size), "clusters\n")
    cat("\nCluster centers:\n")
    print(x$centers, ...)
    cat("\nMemberships:\n")
    print(x$membership, ...)
    cat("\nClosest hard clustering:\n")
    print(x$cluster, ...)
    cat("\nAvailable components:\n")
    print(names(x), ...)
    invisible(x)
}
