# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## 'simApply' returns data.frame with information about design 

setMethod("simApply", 
    signature(x = "data.frame", design = "BasicVector", fun = "function"), 
    function(x, design, fun, ...) {
        s <- stratify(x, design)
        simApply(x, s, fun, ...)
    })

# the function to be applied must return a vector 
setMethod("simApply", 
    signature(x = "data.frame", design = "Strata", fun = "function"), 
    function(x, design, fun, ...) {
        xSpl <- lapply(getSplit(design), function(s, x) x[s,], x)
        tmp <- lapply(xSpl, 
            function(x, ...) {
                res <- fun(x, ...)
                if(is.null(res) || is.vector(res)) res
                else stop("'fun' must return a vector")
            }, ...)
        ind <- sapply(tmp, function(x) as.logical(length(x)))
        values <- do.call("rbind", tmp)
        data.frame(getLegend(design)[ind, , drop=FALSE], values)
    })

# the function to be applied may return a vector, matrix or data.frame
# (not be desirable for our purposes, at least not now)
#setMethod("simApply", 
#    signature(x = "data.frame", design = "Strata", fun = "function"), 
#    function(x, design, fun, ...) {
#        xSpl <- lapply(getSplit(design), function(s, x) x[s,], x)
#        tmp <- lapply(xSpl, 
#            function(x, ...) {
#                res <- fun(x, ...)
#                if(is.null(res) || is.vector(res) || 
#                    is.matrix(res) || is.data.frame(res)) res
#                else stop("'fun' must return a vector, matrix or data.frame")
#            }, ...)
#        nrep <- sapply(tmp, 
#            function(x) {
#                if(is.null(x)) 0
#                else if(is.vector(x)) {
#                    if(length(x)) 1 else 0
#                } else nrow(x)
#            })
#        ind <- rep(getNr(design), each=nrep)
#        values <- do.call("rbind", tmp)
#        res <- data.frame(getLegend(design)[ind, , drop=FALSE], values)
#        row.names(res) <- NULL
#        res
#    })
