### * cl_proximity

cl_proximity <-
function(x, description, class = NULL, labels = NULL,
         self = NULL, size = NULL)
{
    ## Similar to as.dist(), in a way.
    ## Currently, as.dist() is not generic, so we cannot provide a
    ## cl_proximity method for it.  Hence, we have our dissimilarities
    ## and ultrametrics extend dist, and we use capitalized names for
    ## the attributes provided for compatibility with dist (Size and
    ## Labels).
    
    if(inherits(x, "dist")) {
        ## Explicitly deal with dist objects.
        ## Useful in particular because cophenetic() returns them.
        out <- x
        if(is.null(size))
            size <- attr(x, "Size")
        if(is.null(labels))
            labels <- attr(x, "Labels")
    }
    else if(inherits(x, "cl_proximity") ||
            !(is.matrix(x) && (nrow(x) == ncol(x))))
        out <- x
    else {
        ## Actually, x should really be a square symmetric matrix.
        ## The "self-proximities" in the main diagonal must be stored
        ## provided there is one non-zero entry.
        self <- diag(x)
        if(all(self == 0)) self <- NULL
        out <- x[row(x) > col(x)]
        if(is.null(labels)) {
            if(!is.null(rownames(x))) 
                labels <- rownames(x)
            else if(!is.null(colnames(x))) 
                labels <- colnames(x)
        }
    }
    if(is.null(size))
        size <- as.integer((sqrt(1 + 8 * length(out)) + 1) / 2)
    attributes(out) <-
        list(Size = size, Labels = labels, description = description,
             self = self)
    class(out) <- unique(c(class, "cl_proximity"))
    out
}

### * names.cl_proximity

names.cl_proximity <-
function(x)
    NULL

### * print.cl_proximity

print.cl_proximity <-
function(x, ...)
{
    description <- attr(x, "description")
    if(length(description) > 0L) {
        ## Could make this generic ...
        kind <- if(inherits(x, "cl_dissimilarity"))
            "Dissimilarities"
        else if(inherits(x, "cl_agreement"))
            "Agreements"
        else
            "Proximities"
        cat(sprintf("%s using %s", kind, description), ":\n", sep = "")
    }
    m <- format(as.matrix(x))
    if(is.null(self <- attr(x, "self")))
        m[row(m) <= col(m)] <- ""
    else
        m[row(m) < col(m)] <- ""
    print(if(is.null(self)) m[-1, -attr(x, "Size")] else m,
          quote = FALSE, right = TRUE, ...)
    invisible(x)
}

### * as.matrix.cl_proximity

as.matrix.cl_proximity <-
function(x, ...)
{
    size <- attr(x, "Size")
    m <- matrix(0, size, size)
    m[row(m) > col(m)] <- x
    m <- m + t(m)
    if(!is.null(self <- attr(x, "self"))) {
        diag(m) <- self
    }
    ## <NOTE>
    ## stats:::as.matrix.dist() provides default dimnames
    ## (seq_len(size)) if no labels are available.
    ## We used to do this too, but ...
    if(!is.null(labels <- attr(x, "Labels")))
        dimnames(m) <- list(labels, labels)
    ## </NOTE>
    m
}

### * [.cl_proximity

"[.cl_proximity" <-
function(x, i, j)
{
    ## Subscripting proximity objects.
    ## Basically matrix-like, but proximity objects are always
    ## "matrices", hence no 'drop' argument.
    ## For double-index subscripting, if i and j are identical,
    ## structure and class are preserved.  Otherwise, a cross-proximity
    ## object is returned (and methods for classes inheriting from
    ## proximity need to readjust the class info as needed).
    ## For single-index subscripting, no attempty is currently made at
    ## preserving structure and class where possible.  (We might also
    ## change this to select objects, i.e., the same rows and columns.)
    size <- attr(x, "Size")
    if(missing(j)) {
        if(missing(i))
            return(x)
        else
            j <- seq_len(size)
    }
    if(missing(i))
        i <- seq_len(size)
    description <- attr(x, "description")
    ## RG's graph:::[.dist avoids as.matrix() in noting that for dist
    ## objects, entry (i,j) is at n(i-1) - i(i-1)/2 + j - i (in the
    ## veclh dist representation).  We could do something similar, but
    ## note that not all proximities have zero diagonals (i.e., NULL
    ## "self" attributes).
    y <- as.matrix(x)[i, j, drop = FALSE]
    if(identical(i, j)) {
        ## Testing using identical() is rather defensive ...
        return(cl_proximity(y, description = description,
                            class = class(x)))
    }
    cl_cross_proximity(y, description = description)
}

### * cl_cross_proximity

cl_cross_proximity <-
function(x, description = NULL, class = NULL)
{
    attr(x, "description") <- description
    class(x) <- c(class, "cl_cross_proximity")
    x
}

### * print.cl_cross_proximity

print.cl_cross_proximity <-
function(x, ...)
{
    description <- attr(x, "description")
    if(length(description) > 0L) {
        ## Could make this generic ...
        kind <- if(inherits(x, "cl_cross_dissimilarity"))
            "Cross-dissimilarities"
        else if(inherits(x, "cl_cross_agreement"))
            "Cross-agreements"
        else
            "Cross-proximities"
        cat(sprintf("%s using %s", kind, description), ":\n", sep = "")
    }
    print(matrix(as.vector(x), nrow = nrow(x), dimnames = dimnames(x)),
          ...)
    invisible(x)
}

### ** print_description_prefix

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
