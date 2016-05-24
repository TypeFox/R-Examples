
# knn.r
#
#   implements generic k-nearest neighbors, i.e. for arbitrary distance 
#   measures, in a way that is compatible with "knn" in package class. 
#
# ceeboo 2005

gknn <- function(x, y, k=1, l=0, break.ties=TRUE, use.all=TRUE, prob=FALSE) {
    if (!is.matrix(x))
       stop(paste(sQuote("x"),"not a matrix"))
    if (!is.factor(y))
       stop(paste(sQuote("y"),"not a factor"))
    if (length(y) != dim(x)[2])
       stop(paste(sQuote("x"),"and",sQuote("y"),"non-conformable"))
    storage.mode(x) <- "double"
    storage.mode(y) <- storage.mode(k) <- storage.mode(l) <- "integer"
    storage.mode(break.ties) <- storage.mode(use.all) <- storage.mode(prob) <- "logical"
    #
    y <- .Call(R_gknn, x, y, k, l, break.ties, use.all, prob)
    y
}

### the end
