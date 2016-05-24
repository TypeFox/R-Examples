## Misc. utility functions for bdpopt

## Extract the rows or columns of a matrix as lists
mat.rows <- function(X) lapply(1:nrow(X), function(i) X[i,])
mat.cols <- function(X) lapply(1:ncol(X), function(i) X[,i])
    
## Map a function over the rows or columns of a matrix
map.rows <- function(f, X) vapply(mat.rows(X), f, 0)
map.cols <- function(f, X) vapply(mat.cols(X), f, 0)

## Create the cartesian product of a list (length >= 1) of numeric vectors.
## The first element changes fastest with the implementation below.
cart.prod <- function(vectors) {
    if (length(vectors) == 1)
        lapply(vectors[[1]], identity)
    else {
        first <- vectors[[1]]
        rest <- cart.prod(vectors[-1])
        do.call(c, lapply(rest, function(r) lapply(first, function(x) c(x, r))))
    }
}

## This alternative version using expand.grid is somewhat slower than the recursive approach above
cart.prod2 <- function(vectors) {
    X <- as.matrix(expand.grid(vectors, KEEP.OUT.ATTRS = FALSE))
    colnames(X) <- NULL
    mat.rows(X)  
}

## Check for elementwise all.equality (i.e., equality within a small tolerance)
elementwise.all.equal <- Vectorize(function(x, y) {isTRUE(all.equal(x, y))})

## Compute the trace of a matrix
tr <- function(M) sum(diag(M))

## Convert a vector of names, some referring to array elements, e.g., c("a", "b[1,1]", ..., "b[3,2]") into a list of named dimensions corresponding to the given vector in order.
## So, the output of the example above would be list(a = c(1), b = c(3,2)).
## It is assumed that the array element with the largest index in all dimensions comes last in order.
named.dims.from.elements <- function(elements) {
    named.dims <- list()

    for (s in elements) {
        split <- strsplit(s, "[", fixed = TRUE)[[1]]
        
        if (length(split) == 1)
            named.dims[[split]] <- c(1)
        else {
            indices <- substring(split[2], 1, nchar(split[2]) - 1) ## remove the ] at the end
            named.dims[[split[1]]] <- as.numeric(strsplit(indices, ",", fixed = TRUE)[[1]])
        }
    }

    named.dims
}

## Construct a list of named arrays from a list of named dimension vectors
## and a vector of numbers of total length equal to the sum of dimension vector products.
construct.named.arrays <- function(named.dims, v) {
    named.arrays <- list()
    i <- 1

    for (name in names(named.dims)) {
        dims <- named.dims[[name]]
        prod.dims <- prod(dims)
        named.arrays[[name]] <- array(v[i:(i + prod.dims - 1)], dim = dims)
        i <- i + prod.dims
    }

    named.arrays
}
