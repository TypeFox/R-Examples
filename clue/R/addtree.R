### * ls_fit_addtree

ls_fit_addtree <-
function(x, method = c("SUMT", "IP", "IR"), weights = 1,
         control = list())
{
    if(!inherits(x, "dist"))
        x <- as.dist(x)

    ## Catch some special cases right away.
    if(attr(x, "Size") <= 3L)
        return(as.cl_addtree(x))
    if(.non_additivity(x, max = TRUE) == 0)
        return(as.cl_addtree(x))

    ## Handle argument 'weights'.
    ## This is somewhat tricky ...
    if(is.matrix(weights)) {
        weights <- as.dist(weights)
        if(length(weights) != length(x))
            stop("Argument 'weights' must be compatible with 'x'.")
    }
    else
        weights <- rep(weights, length.out = length(x))
    if(any(weights < 0))
        stop("Argument 'weights' has negative elements.")
    if(!any(weights > 0))
        stop("Argument 'weights' has no positive elements.")

    method <- match.arg(method)
    switch(method,
           SUMT = .ls_fit_addtree_by_SUMT(x, weights, control),
           IP = {
               .ls_fit_addtree_by_iterative_projection(x, weights,
                                                       control)
           },
           IR = {
               .ls_fit_addtree_by_iterative_reduction(x, weights,
                                                      control)
           })
}

### ** .ls_fit_addtree_by_SUMT

.ls_fit_addtree_by_SUMT <-
function(x, weights = 1, control = list())
{
    ## Control parameters:
    ## gradient,
    gradient <- control$gradient
    if(is.null(gradient))
        gradient <- TRUE
    ## nruns,
    nruns <- control$nruns
    ## start,
    start <- control$start

    ## Handle start values and number of runs.
    if(!is.null(start)) {
        if(!is.list(start)) {
            ## Be nice to users.
            start <- list(start)
        }
    } else if(is.null(nruns)) {
        ## Use nruns only if start is not given.
        nruns <- 1L
    }

    w <- weights / sum(weights)
    n <- attr(x, "Size")
    labels <- attr(x, "Labels")    

    ## Handle missing values in x along the lines of de Soete (1984):
    ## set the corresponding weights to 0, and impute by the weighted
    ## mean.
    ind <- which(is.na(x))
    if(any(ind)) {
        w[ind] <- 0
        x[ind] <- weighted.mean(x, w, na.rm = TRUE)
    }

    L <- function(d) sum(w * (d - x) ^ 2)
    P <- .make_penalty_function_addtree(n)
    if(gradient) {
        grad_L <- function(d) 2 * w * (d - x)
        grad_P <- .make_penalty_gradient_addtree(n)
    }
    else {
        grad_L <- grad_P <- NULL
    }

    if(is.null(start)) {
        ## Initialize by "random shaking".  Use sd() for simplicity.
        start <- replicate(nruns,  
                           x + rnorm(length(x), sd = sd(x) / sqrt(3)),
                           simplify = FALSE)
    }

    ## And now ...
    d <- sumt(start, L, P, grad_L, grad_P,
              method = control$method, eps = control$eps,
              q = control$q, verbose = control$verbose,
              control = as.list(control$control))$x

    ## Round to enforce additivity, and hope for the best ...
    .cl_addtree_from_addtree_approximation(d, n, labels)    

}

.make_penalty_function_addtree <-
function(n)
    function(d) {
        (.non_additivity(.symmetric_matrix_from_veclh(d, n))
         + sum(pmin(d, 0) ^ 2))
    }

.make_penalty_gradient_addtree <-
function(n)
    function(d) {
        gr <- matrix(.C(C_deviation_from_additivity_gradient,
                        as.double(.symmetric_matrix_from_veclh(d, n)),
                        as.integer(n),
                        gr = double(n * n))$gr,
                     n, n)
        gr[row(gr) > col(gr)] + 2 * sum(pmin(d, 0))
    }
    

### ** .ls_fit_addtree_by_iterative_projection

## <NOTE>
## Functions
##   .ls_fit_addtree_by_iterative_projection()
##   .ls_fit_addtree_by_iterative_reduction()
## are really identical apart from the name of the C routine they call.
## (But will this necessarily always be the case in the future?)
## Merge maybe ...
## </NOTE>

.ls_fit_addtree_by_iterative_projection <-
function(x, weights = 1, control = list())
{
    if(any(diff(weights)))
        warning("Non-identical weights currently not supported.")

    labels <- attr(x, "Labels")
    x <- as.matrix(x)
    n <- nrow(x)

    ## Control parameters:
    ## maxiter,
    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 10000L
    ## nruns,
    nruns <- control$nruns
    ## order,
    order <- control$order
    ## tol,
    tol <- control$tol
    if(is.null(tol))
        tol <- 1e-8
    ## verbose.
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")

    ## Handle order and nruns.
    if(!is.null(order)) {
        if(!is.list(order))
            order <- as.list(order)
        if(!all(sapply(order,
                       function(o) all(sort(o) == seq_len(n)))))
            stop("All given orders must be valid permutations.")
    }
    else {
        if(is.null(nruns))
            nruns <- 1L
        order <- replicate(nruns, sample(n), simplify = FALSE)
    }

    ind <- lower.tri(x)
    L <- function(d) sum(weights * (x - d)[ind] ^ 2)

    d_opt <- NULL
    v_opt <- Inf
    for(run in seq_along(order)) {
        if(verbose)
            message(gettextf("Iterative projection run: %d", run))
        d <- .C(C_ls_fit_addtree_by_iterative_projection,
                as.double(x),
                as.integer(n),
                as.integer(order[[run]] - 1L),
                as.integer(maxiter),
                iter = integer(1L),
                as.double(tol),
                as.logical(verbose))[[1L]]
        v <- L(d)
        if(v < v_opt) {
            v_opt <- v
            d_opt <- d
        }
    }

    d <- matrix(d_opt, n)
    dimnames(d) <- list(labels, labels)
    .cl_addtree_from_addtree_approximation(as.dist(d))
}

### ** .ls_fit_addtree_by_iterative_reduction

.ls_fit_addtree_by_iterative_reduction <-
function(x, weights = 1, control = list())
{
    if(any(diff(weights)))
        warning("Non-identical weights currently not supported.")

    labels <- attr(x, "Labels")
    x <- as.matrix(x)
    n <- nrow(x)

    ## Control parameters:
    ## maxiter,
    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 10000L
    ## nruns,
    nruns <- control$nruns
    ## order,
    order <- control$order
    ## tol,
    tol <- control$tol
    if(is.null(tol))
        tol <- 1e-8
    ## verbose.
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")

    ## Handle order and nruns.
    if(!is.null(order)) {
        if(!is.list(order))
            order <- as.list(order)
        if(!all(sapply(order,
                       function(o) all(sort(o) == seq_len(n)))))
            stop("All given orders must be valid permutations.")
    }
    else {
        if(is.null(nruns))
            nruns <- 1L
        order <- replicate(nruns, sample(n), simplify = FALSE)
    }

    ind <- lower.tri(x)
    L <- function(d) sum(weights * (x - d)[ind] ^ 2)

    d_opt <- NULL
    v_opt <- Inf
    for(run in seq_along(order)) {
        if(verbose)
            message(gettextf("Iterative reduction run: %d", run))
        d <- .C(C_ls_fit_addtree_by_iterative_reduction,
                as.double(x),
                as.integer(n),
                as.integer(order[[run]] - 1L),
                as.integer(maxiter),
                iter = integer(1L),
                as.double(tol),
                as.logical(verbose))[[1L]]
        v <- L(d)
        if(v < v_opt) {
            v_opt <- v
            d_opt <- d
        }
    }

    d <- matrix(d_opt, n)
    dimnames(d) <- list(labels, labels)
    .cl_addtree_from_addtree_approximation(as.dist(d))
}

### * .non_additivity

.non_additivity <-
function(x, max = FALSE)
{
    if(!is.matrix(x))
        x <- .symmetric_matrix_from_veclh(x)
    .C(C_deviation_from_additivity,
       as.double(x),
       as.integer(nrow(x)),
       fn = double(1L),
       as.logical(max))$fn
}

### * ls_fit_centroid

ls_fit_centroid <-
function(x)
{
    ## Fit a centroid additive tree distance along the lines of Carroll
    ## & Pruzansky (1980).  In fact, solving
    ##
    ##    \sum_{i,j: i \ne j} (\delta_{ij} - (g_i + g_j)) ^ 2 => min_g
    ##
    ## gives \sum_{j: j \ne i} (g_i + g_j - \delta_{ij}) = 0, or (also
    ## in Barthemely & Guenoche)
    ##
    ##    (n - 2) g_i + \sum_j g_j = \sum_{j: j \ne i} \delta_{ij}
    ##
    ## which after summing over all i and some manipulations eventually
    ## gives
    ##
    ##    g_i = \frac{1}{n-2} (v_i - m),
    ##
    ##    v_i = \sum_{j: j \ne i} \delta_{ij}
    ##    s   = \frac{1}{2(n-1)} \sum_{i,j: j \ne i} \delta_{ij}

    n <- attr(x, "Size")

    if(n <= 2L)
        return(as.cl_addtree(0 * x))

    x <- as.matrix(x)
    g <- rowSums(x) / (n - 2) - sum(x) / (2 * (n - 1) * (n - 2))
    as.cl_addtree(as.dist(.make_centroid_matrix(g)))
}

.make_centroid_matrix <-
function(g)
{
    y <- outer(g, g, "+")
    diag(y) <- 0
    y
}

### * as.cl_addtree

as.cl_addtree <-
function(x)
    UseMethod("as.cl_addtree")

as.cl_addtree.default <-
function(x)    
{
    if(inherits(x, "cl_addtree"))
        x
    else if(is.atomic(x) || inherits(x, "cl_ultrametric"))
        .cl_addtree_from_veclh(x)
    else if(is.matrix(x)) {
        ## Should actually check whether the matrix is symmetric, >= 0
        ## and satisfies the 4-point conditions ...
        .cl_addtree_from_veclh(as.dist(x))
    }
    else if(is.cl_dendrogram(x))
        .cl_addtree_from_veclh(cl_ultrametric(x))
    else
        stop("Cannot coerce to 'cl_addtree'.")
}

as.cl_addtree.phylo <-
function(x)
    .cl_addtree_from_veclh(as.dist(cophenetic(x)))
## Phylogenetic trees with edge/branch lengths yield additive tree
## dissimilarities.

### * .cl_addtree_from_veclh

.cl_addtree_from_veclh <-
function(x, size = NULL, labels = NULL)
{
    cl_proximity(x, "Additive tree distances",
                 labels = labels, size = size,
                 class = c("cl_addtree", "cl_dissimilarity",
                 "cl_proximity", "dist"))
}

### * .cl_addtree_from_addtree_approximation

.cl_addtree_from_addtree_approximation <-
function(x, size = NULL, labels = NULL)
{
    ## Turn x into an addtree after possibly rounding to non-additivity
    ## significance (note that this is not guaranteed to work ...).
    mnum <- .non_additivity(x, max = TRUE)
    x <- round(x, floor(abs(log10(mnum))))
    .cl_addtree_from_veclh(x, size = size, labels = labels)
}

### * .decompose_addtree

.decompose_addtree <-
function(x, const = NULL)
{
    ## Decompose an addtree into an ultrametric and a centroid
    ## distance.

    ## If 'const' is not given, we take the root as half way between the
    ## diameter of the addtree, and choose a minimal constant to ensure
    ## non-negativity (but not positivity) of the ultrametric.

    ## As this is all slightly dubious and it is not quite clear how
    ## much positivity we want in the ultrametric of the decomposition,
    ## we keep this hidden.  For plotting addtrees, the choice of the
    ## constant does not seem to matter.

    x <- as.matrix(x)
    n <- nrow(x)

    ## Determine diameter.
    ind <- which.max(x) - 1
    u <- ind %% n + 1
    v <- ind %/% n + 1

    if(!is.null(const))
        g <- pmax(x[u, ], x[v, ]) - const
    else {
        g <- pmax(x[u, ], x[v, ]) - x[u, v] / 2
        u <- x - .make_centroid_matrix(g)
        k <- - min(u)
        g <- g - k / 2
    }
    u <- x - .make_centroid_matrix(g)

    names(g) <- rownames(x)

    ## Ensure a valid ultrametric.
    d <- .ultrametrify(as.dist(u))
    u <- .cl_ultrametric_from_veclh(d, nrow(x), rownames(x))

    ## Note that we return the centroid distances to the root, and not
    ## between the objects (as.dist(.make_centroid_matrix(g))) ...
    list(Ultrametric = as.cl_ultrametric(u), Centroid = g)
}

### * plot.cl_addtree

plot.cl_addtree <-
function(x, ...)
{
    ## Construct a dendrogram-style representation of the addtree with
    ## the root half way between the diameter, and plot.
    y <- .decompose_addtree(x, max(x))
    u <- y$Ultrametric
    g <- y$Centroid
    ## We halve the scale of the ultrametric, and add the maximal g from
    ## the centroid.
    h <- hclust(as.dist(u / 2), "single")
    h$height <- h$height + max(g)
    d <- as.dendrogram(h)
    ## Now modify the heights of the leaves so that the objects giving
    ## the diameter of the addtree end up with height zero.
    g <- max(g) - g
    names(g) <- labels(g)
    d <- dendrapply(d, function(n) {
        if(!is.leaf(n)) return(n)
        attr(n, "height") <- g[attr(n, "label")]
        n
    })
    ## And finally plot
    plot(d, ...)
}


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
