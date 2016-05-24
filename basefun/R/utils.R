
### Box product of two matrices
.boxprod <- function(x, y) {
    ret <- x[, rep(1:ncol(x), ncol(y)), drop = FALSE] * 
           y[, rep(1:ncol(y), rep(ncol(x), ncol(y))), drop = FALSE]
    colnames(ret) <- paste(colnames(x)[rep(1:ncol(x), ncol(y))], 
                           colnames(y)[rep(1:ncol(y), rep(ncol(x), ncol(y)))], sep = ":")
    ret
}

### box product of design matices
.box <- function(...) {
    args <- list(...)
    ret <- args[[1]]
    if (length(args) == 1) return(ret)
    for (i in 2:length(args))
        ret <- .boxprod(ret, args[[i]])
    ret
}

### Box "product" of two character matrices
.boxprod_char <- function(x, y) {
    matrix(paste(x[rep(1:nrow(x), length = nrow(y)), 
                   rep(1:ncol(x), ncol(y)), drop = FALSE], 
                 y[rep(1:nrow(y), length = nrow(x)), 
                   rep(1:ncol(y), rep(ncol(x), ncol(y))), drop = FALSE], 
                 sep = ":"), nrow = nrow(x))
}

### box "product" of character matrices
.box_char <- function(...) {
    args <- list(...)
    ret <- args[[1]]
    if (length(args) == 1) return(ret)
    for (i in 2:length(args))
        ret <- .boxprod_char(ret, args[[i]])
    ret
}

### linear functions for constraints corresponding to box products
.box_ui_ci <- function(...) {
    args <- list(...)
    ret <- args[[1]]
    if (length(args) == 1) return(ret)
    for (i in 2:length(args)) {
        nci <- ncol(args[[i]]$ui)
        ncr <- ncol(ret$ui)
        ret$ui <- rBind(kronecker(Diagonal(nci), ret$ui),
                        kronecker(args[[i]]$ui, Diagonal(ncr)))
        ret$ci <- c(as(kronecker(Diagonal(nci), 
                                 matrix(ret$ci, ncol = 1)) %*% rep(1, nci), 
                       "vector"),
                    as(kronecker(matrix(args[[i]]$ci, ncol = 1), 
                                 Diagonal(ncr)) %*% rep(1, ncr), 
                       "vector"))
    }
    ret
}

.const_array <- function(dim, subdim, x) {
    nd <- names(dim)
    stopifnot(all(subdim %in% nd))
    stopifnot(length(x) == prod(dim[subdim]))
    nd2 <- c(nd[nd %in% subdim], nd[!(nd %in% subdim)])
    
    ret <- array(x, dim = dim[nd2])
    aperm(ret, perm = match(nd, nd2))

}

.deriv <- function(vn, deriv) {

    ### applies to all basis functions
    ### deriv = -1 => 0
    ### deriv = 0 => f(x)
    ### deriv = 1 => f'(x)
    ### deriv = 2 => f''(x)
    ### ...

    stopifnot(length(deriv) == 1)
    if (deriv == 0L) return(deriv)
    stopifnot(!is.null(names(deriv)))
    if (names(deriv) %in% vn)
        return(deriv)
    return(-1L)
}

L2B <- function(order) {
    stopifnot(order > 0)
    .Call("L2B", as.integer(order))
}

B2L <- function(order) {
    stopifnot(order > 0)
    .Call("B2L", as.integer(order))
}
