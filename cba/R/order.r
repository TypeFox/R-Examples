
# wrapper to the optimal leaf ordering algorithm
#
# ceeboo 2005

order.optimal <- function(dist, merge) {
    if (!inherits(dist,"dist"))
       stop(paste(sQuote("dist"),"not of class dist"))
    if (!is.matrix(merge))
       stop(paste(sQuote("merge"),"not a matrix"))
    if (length(dim(merge)) != 2)
       stop(paste(sQuote("merge"),"invalid"))
    if (dim(merge)[1] != attr(dist,"Size")-1)
       stop(paste(sQuote("dist"),"and",sQuote("merge"),"do not conform"))
    if (!is.double(dist))
       storage.mode(dist) <- "double"
    storage.mode(merge) <- "integer"
    obj <- .Call(R_order_optimal, dist, merge)
    names(obj) <- c("merge","order","length")
    names(obj$order) <- attr(dist,"Labels")
    obj
}

# wrapper to computing the lenght of the order
# under a distance matrix, e.g. a tour where the
# leg between the first and last city is omitted.
# that this is a (Hamilton) path.
#
# note that this corresponds to the sum of distances 
# along the first off diagonal of the ordered distance
# matrix.
# 

order.length <- function(dist, order) {
    if (!inherits(dist,"dist"))
       stop(paste(sQuote("dist"),"not of class dist"))
    if (missing(order))
       order <- 1:attr(dist, "Size")
    else {
       if (length(order) != attr(dist,"Size"))
          stop(paste(sQuote("order"),"invalid lenght"))
    }
    if (!is.double(dist))
       storage.mode(dist) <- "double"
    if (!is.integer(order))
       storage.mode(order) <- "integer"
    x <- .Call(R_order_length, dist, order)
    x
}

# wrapper to greedy ordering inspired by F. Murtagh
# actually a hierarchical cluster algorithm.

order.greedy <- function(dist) {
    if (!inherits(dist, "dist"))
       stop(paste(sQuote("dist"),"not of class dist"))
    if (!is.double(dist))
       storage.mode(dist) <- "double"
    obj <- .Call(R_order_greedy, dist)
    names(obj) <- c("merge", "order", "height");
    obj
}

###
