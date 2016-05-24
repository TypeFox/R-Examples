### * Matrix/vector utilities

### * .dist_from_vector
        
.dist_from_vector <-
function(x, n = NULL, labels = NULL)
{
    ## This might be useful as as.dist.vector, perhaps without the extra
    ## argument n then which we only have for minimal performance gains.
    if(is.null(n))
        n <- as.integer((sqrt(1 + 8 * length(x)) + 1) / 2)
    attr(x, "Size") <- n
    if(!is.null(labels))
        attr(x, "Labels") <- labels
    class(x) <- "dist"
    x
}

### ** .one_entry_per_column

.one_entry_per_column <-
function(x, j)
{
    ## For a matrix x and a vector of column indices j_1, ..., j_n where
    ## n is the number of rows of x, get x[1,j_1], ..., x[n,j_n].

    ## <NOTE>
    ## This used to have
    ##   if(!is.matrix(x))
    ##     stop("Argument 'x' must be a matrix.")
    ## but that will fail for sparse matrix classes.
    ## So let us hope for the best ...
    ## </NOTE>

    x[cbind(seq_len(nrow(x)), j)]
}

".one_entry_per_column<-" <-
function(x, j, value)
{
    ## <NOTE>
    ## This used to have
    ##   if(!is.matrix(x))
    ##     stop("Argument 'x' must be a matrix.")
    ## but that will fail for sparse matrix classes.
    ## So let us hope for the best ...
    ## </NOTE>

    x[cbind(seq_len(nrow(x)), j)] <- value
    x
}

### * .symmetric_matrix_from_veclh

.symmetric_matrix_from_veclh <-
function(x, n = NULL)
{
    ## In essence the same as as.matrix.dist, but without handling the
    ## additional attributes that dist objects might have.
    if(is.null(n))
        n <- as.integer((sqrt(1 + 8 * length(x)) + 1) / 2)
    M <- matrix(0, n, n)
    M[row(M) > col(M)] <- x
    M + t(M)
}

### * .weighted_mean_of_object_dissimilarities

.weighted_mean_of_object_dissimilarities <-
function(x, w = NULL)
{
    w <- if(is.null(w)) {
        rep.int(1, length(x))
    } else {
        rep(w, length.out = length(x))
    }
    ## (Need the latter because we want w / sum(w) ...)
    dissimilarities <- lapply(x, cl_object_dissimilarities)
    m <- rowSums(mapply("*", dissimilarities, w / sum(w)))
    labels <- attr(dissimilarities[[1L]], "Labels")
    .dist_from_vector(m, labels = labels)
}

### ** .weighted_sum_of_matrices

.weighted_sum_of_matrices <-
function(x, w = NULL, nr = NULL)
{
    ## Quite often we need to compute weighted sums \sum_b w_b X_b of
    ## conforming matrices \{ X_b \}.  If x is a list containing the
    ## matrices and w the vector of weights, it seems that one
    ## reasonably efficient way of doing this is the following.
    if(is.null(w)) w <- rep.int(1, length(x))
    if(is.null(nr)) nr <- NROW(x[[1L]])
    matrix(rowSums(mapply("*", x, w)), nr)
}

### ** .weighted_sum_of_vectors

.weighted_sum_of_vectors <-
function(x, w = NULL)
{
    ## See above.
    if(is.null(w)) w <- rep.int(1, length(x))
    rowSums(mapply("*", x, w))
}

### * Containers

## Creator.

.make_container <-
function(x, classes, properties = NULL)
{
    out <- list(.Data = x, .Meta = properties)
    class(out) <- unique(classes)
    out
}
          
## Getters.

.get_representation <-
function(x)
    x$.Data
.get_properties <-
function(x)
    x$.Meta
.get_property <-
function(x, which)
    x$.Meta[[which]]
.has_property <-
function(x, which)
    which %in% names(x$.Meta)
.get_property_from_object_or_representation <-
function(x, which, getter)
{
    if(.has_property(x, which))
        .get_property(x, which)
    else {
        if(missing(getter)) getter <- get(which)
        getter(.get_representation(x))
    }
}

## Methods (sort of).

.print_container <-
function(x, cls, ...)
{
    writeLines(gettextf("An object of virtual class '%s', with representation:\n",
                        cls))
    print(.get_representation(x), ...)
    invisible(x)
}
    
### * Others

weighted_median <-
function(x, w = 1, na.rm = FALSE)
{
    w <- rep(w, length.out = length(x))
    if(na.rm && any(ind <- is.na(x))) {
        x <- x[!ind]
        w <- w[!ind]
    }
    if(any(is.na(x)) || !length(x))
        return(NA)
    w <- w / sum(w)    
    ind <- order(x)
    x <- x[ind]
    w <- w[ind]
    x[which.min(x * (cumsum(w) - 0.5) - cumsum(w * x))]
}


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
