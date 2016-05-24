### * cl_ultrametric

cl_ultrametric <-
function(x, size = NULL, labels = NULL) 
{
    if(inherits(x, "cl_hierarchy")) {
        ## <NOTE>
        ## Strictly, not every hierarchy corresponds to an ultrametric.
        ## </NOTE>
        return(cl_ultrametric(.get_representation(x),
                              size = size, labels = labels))
    }
    else if(!inherits(x, "cl_ultrametric")) {
        ## Try using cophenetic().
        ## This starts by coercing to hclust, which has methods for all
        ## currently supported hierarchical classification methods.
        ## To support others, either provide as.hclust methods for
        ## these, or make cl_ultrametric() generic and add methods.
        ## Or use the fact that in R >= 2.1.0, stats::cophenetic() is
        ## generic.
        out <- cophenetic(x)
    }
    else {
        out <- x
        if(is.null(labels))
            labels <- attr(x, "Labels")
    }
    .cl_ultrametric_from_veclh(out, labels = labels, size = size)
}

.cl_ultrametric_from_veclh <-
function(x, size = NULL, labels = NULL, meta = NULL)
{
    if(.non_ultrametricity(x) > 0)
        stop("Not a valid ultrametric.")
    u <- cl_proximity(x, "Ultrametric distances",
                      labels = labels, size = size,
                      class = c("cl_ultrametric", "cl_dissimilarity",
                      "cl_proximity", "dist"))
    if(!is.null(meta))
        attr(u, "meta") <- meta
    u
}

### * as.cl_ultrametric

as.cl_ultrametric <-
function(x)
    UseMethod("as.cl_ultrametric")
as.cl_ultrametric.default <-
function(x)
{
    if(inherits(x, "cl_ultrametric"))
        x
    else if(is.atomic(x))
        .cl_ultrametric_from_veclh(x)
    else
        cl_ultrametric(x)
}
as.cl_ultrametric.matrix <-
function(x)
    .cl_ultrametric_from_veclh(x[row(x) > col(x)],
                               labels = rownames(x))

### * as.dendrogram.cl_ultrametric

as.dendrogram.cl_ultrametric <-
function(object, ...)
    as.dendrogram(as.hclust(object), ...)

### * as.hclust.cl_ultrametric

as.hclust.cl_ultrametric <-
function(x, ...)
{
    ## Hierarchical clustering with single linkage gives the minimal
    ## ultrametric dominated by a dissimilarity, see e.g. Bock (1974,
    ## Theorem 39.2).  Hence, hclust(method = "single") on an
    ## ultrametric gives the hclust representation of the associated
    ## dendrogram.
    hclust(x, "single")
}

### * cophenetic.cl_ultrametric

cophenetic.cl_ultrametric <-
function(x)
    as.dist(x)

### * plot.cl_ultrametric

plot.cl_ultrametric <-
function(x, ...)
    plot(as.dendrogram(x), ...)

### * ls_fit_ultrametric

ls_fit_ultrametric <-
function(x, method = c("SUMT", "IP", "IR"), weights = 1,
         control = list())
{
    if(inherits(x, "cl_ultrametric")) {
        return(.cl_ultrametric_with_meta_added(x, list(objval = 0)))
    } else if(is.cl_ensemble(x) || is.list(x)) {
        ## Might be given a list/ensemble of object dissimilarities.
        ## In this case, compute the suitably weighted average and
        ## proceed.
        if(length(x) == 0L)
            stop("Given ensemble contains no dissimilarities.")
        ## Let's be nice as usual ...
        ind <- !sapply(x, .has_object_dissimilarities)
        if(any(ind))
            x[ind] <- lapply(x[ind], as.dist)
        x <- .weighted_mean_of_object_dissimilarities(x,
                                                      control$weights)
    } else if(!inherits(x, "dist"))
        x <- as.dist(x)

    ## Catch some special cases right away.
    if(attr(x, "Size") <= 2L)
        return(.cl_ultrametric_with_meta_added(as.cl_ultrametric(x),
                                               list(objval = 0)))
    if(.non_ultrametricity(x, max = TRUE) == 0)
        return(.cl_ultrametric_with_meta_added(as.cl_ultrametric(x),
                                               list(objval = 0)))

    ## Handle weights.
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
           SUMT = .ls_fit_ultrametric_by_SUMT(x, weights, control),
           IP = {
               .ls_fit_ultrametric_by_iterative_projection(x, weights,
                                                           control)
           },
           IR = {
               .ls_fit_ultrametric_by_iterative_reduction(x, weights,
                                                          control)
           })
}

### ** .ls_fit_ultrametric_by_SUMT
           
.ls_fit_ultrametric_by_SUMT <-
function(x, weights = 1, control = list())
{
    ## Fit an ultrametric to a dissimilarity by minimizing euclidean
    ## dissimilarity subject to the ultrametric constraint, using the
    ## sequential algorithm of de Soete (1984) with a slight change: we
    ## try to ensure that what we obtain satisfies the constraints
    ## "exactly" rather than approximately.  We (currently?) do that via
    ## rounding ...

    ## <NOTE>
    ## This fits and hence returns an ultrametric, *not* the hierarchy
    ## corresponding to the ultrametric.
    ## </NOTE>

    w <- weights / sum(weights)

    ## Control parameters:
    ## nruns,
    nruns <- control$nruns
    ## start.
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

    ## If x is an ultrametric, or satisfies the ultrametricity
    ## constraints, return it.
    if(inherits(x, "cl_ultrametric")
       || (.non_ultrametricity(x, max = TRUE) == 0))
        return(.cl_ultrametric_with_meta_added(as.cl_ultrametric(x),
                                               list(objval = 0)))

    ## For the time being, use a simple minimizer.

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
    
    ## We follow de Soete's notation, and use the veclh's (vector of
    ## lower half, in S the same as x[lower.tri(x)]) of the respective
    ## proximity objects.

    L <- function(d) sum(w * (x - d) ^ 2)
    P <- .make_penalty_function_ultrametric(n)
    grad_L <- function(d) 2 * w * (d - x)
    grad_P <- .make_penalty_gradient_ultrametric(n)

    if(is.null(start)) {
        ## Initialize by "random shaking".  Use sd() for simplicity.
        start <- replicate(nruns,  
                           x + rnorm(length(x), sd = sd(x) / sqrt(3)),
                           simplify = FALSE)
    }
    
    ## And now ...
    out <- sumt(start, L, P, grad_L, grad_P,
                method = control$method, eps = control$eps,
                q = control$q, verbose = control$verbose,
                control = as.list(control$control))

    d <- .ultrametrify(out$x)
    meta <- list(objval = L(d))

    .cl_ultrametric_from_veclh(d, n, labels, meta)
}

.make_penalty_function_ultrametric <-
function(n)
    function(d) {
        ## Smooth penalty function measuring the extent of violation of
        ## the ultrametricity constraint.  Also ensure nonnegativity ...
        (.non_ultrametricity(.symmetric_matrix_from_veclh(d, n))
         + sum(pmin(d, 0) ^ 2))
    }

.make_penalty_gradient_ultrametric <-
function(n)
    function(d) {
        gr <- matrix(.C(C_deviation_from_ultrametricity_gradient,
                        as.double(.symmetric_matrix_from_veclh(d, n)),
                        as.integer(n),
                        gr = double(n * n))$gr,
                     n, n)
        gr[row(gr) > col(gr)] + 2 * sum(pmin(d, 0))
    }

### ** .ls_fit_ultrametric_by_iterative_projection

## <NOTE>
## Functions
##   .ls_fit_ultrametric_by_iterative_projection()
##   .ls_fit_ultrametric_by_iterative_reduction()
## are really identical apart from the name of the C routine they call.
## (But will this necessarily always be the case in the future?)
## Merge maybe ...
## </NOTE>

.ls_fit_ultrametric_by_iterative_projection <-
function(x, weights = 1, control = list())
{
    if(any(diff(weights) != 0))
        warning("Non-identical weights currently not supported.")

    labels <- attr(x, "Labels")
    n <- attr(x, "Size")

    x <- as.matrix(x)

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

    ## <NOTE>
    ## Adjust in case support for non-identical weights is added.
    L <- function(d) sum((x - d) ^ 2)
    ## </NOTE>

    d_opt <- NULL
    v_opt <- Inf
    for(run in seq_along(order)) {
        if(verbose)
            message(gettextf("Iterative projection run: %d", run))
        d <- .C(C_ls_fit_ultrametric_by_iterative_projection,
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
        
    d <- .ultrametrify(as.dist(matrix(d_opt, n)))
    meta <- list(objval = L(d))
    
    .cl_ultrametric_from_veclh(d, n, labels, meta)
}

### ** .ls_fit_ultrametric_by_iterative_reduction

.ls_fit_ultrametric_by_iterative_reduction <-
function(x, weights = 1, control = list())
{
    if(any(diff(weights) != 0))
        warning("Non-identical weights currently not supported.")

    labels <- attr(x, "Labels")
    n <- attr(x, "Size")

    x <- as.matrix(x)

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

    ## <NOTE>
    ## Adjust in case support for non-identical weights is added.
    L <- function(d) sum((x - d) ^ 2)
    ## </NOTE>

    d_opt <- NULL
    v_opt <- Inf
    for(run in seq_along(order)) {
        if(verbose)
            message(gettextf("Iterative reduction run: %d", run))
        d <- .C(C_ls_fit_ultrametric_by_iterative_reduction,
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

    d <- .ultrametrify(as.dist(matrix(d_opt, n)))
    meta <- list(objval = L(d))
    
    .cl_ultrametric_from_veclh(d, n, labels, meta)
}

### * Ultrametric Target Fitters.

### ** ls_fit_ultrametric_target

ls_fit_ultrametric_target <-
function(x, y, weights = 1)
{
    fitter <- if(identical(weights, 1)) # Default.
        function(x, w) mean(x)
    else
        function(x, w) weighted.mean(x, w)
    distfun <- function(x, u, w) sqrt(sum(w * (x - u) ^ 2))
    .fit_ultrametric_target(x, y, weights, fitter, distfun)
}

### ** l1_fit_ultrametric_target

l1_fit_ultrametric_target <-
function(x, y, weights = 1)
{
    fitter <- if(identical(weights, 1)) # Default.
        function(x, w) median(x)
    else
        function(x, w) weighted_median(x, w)
    distfun <- function(x, u, w) sum(w * abs(x - u))
    .fit_ultrametric_target(x, y, weights, fitter, distfun)
}

### ** .fit_ultrametric_target    

.fit_ultrametric_target <-
function(x, y, w, fitter, distfun = NULL)
{
    w <- .handle_weights_for_ultrametric_target_fitters(w, x)
    ## The documentation says that x should inherit from dist, so coerce
    ## to this if needed but if not a matrix (as we will coerce back to
    ## a matrix right away).
    if(!inherits(x, "dist") && !is.matrix(x))
        x <- as.dist(x)
    x <- as.matrix(x)
    y <- as.hclust(y)
    n <- length(y$order)
    ilist <- vector("list", n)
    out <- matrix(0, n, n)
    mat <- xlist <- wlist <- vector("list", n - 1L)
    for(i in seq_len(n - 1L)) {
        inds <- y$merge[i, ]
        ids1 <- if(inds[1L] < 0) -inds[1L] else ilist[[inds[1L]]]
        ids2 <- if(inds[2L] < 0) -inds[2L] else ilist[[inds[2L]]]
        ilist[[i]] <- c(ids1, ids2)
        mat[[i]] <- cbind(rep.int(ids1,
                                  rep.int(length(ids2), length(ids1))), 
                          rep.int(ids2, length(ids1)))
        xlist[[i]] <- x[mat[[i]]]
        wlist[[i]] <- w[mat[[i]]]
    }
    values <- pava(xlist, wlist, fitter)
    for(i in seq_len(n - 1L))
        out[mat[[i]]] <- values[i]
    rownames(out) <- y$labels
    u <- as.cl_ultrametric(out + t(out))
    if(!is.null(distfun))
        attr(u, "meta") <-
            list(objval = distfun(as.dist(x), u, as.dist(w)))
    u
}

### ** .handle_weights_for_ultrametric_target_fitters

.handle_weights_for_ultrametric_target_fitters <-
function(weights, x)
{
    ## Handle weights for the ultrametric target fitters.
    ## This is somewhat tricky ...
    if(is.matrix(weights)) {
        if(any(dim(weights) != attr(x, "Size")))
            stop("Argument 'weights' must be compatible with 'x'.")
    }
    else
        weights <-
            as.matrix(.dist_from_vector(rep(weights,
                                            length.out = length(x))))
    if(any(weights < 0))
        stop("Argument 'weights' has negative elements.")
    if(!any(weights > 0))
        stop("Argument 'weights' has no positive elements.")
    weights
}

### l1_fit_ultrametric

l1_fit_ultrametric <-
function(x, method = c("SUMT", "IRIP"), weights = 1,
         control = list())
{
    if(inherits(x, "cl_ultrametric"))
        return(.cl_ultrametric_with_meta_added(x, list(objval = 0)))

    if(!inherits(x, "dist"))
        x <- as.dist(x)

    ## Catch some special cases right away.
    if(attr(x, "Size") <= 2L)
        return(.cl_ultrametric_with_meta_added(as.cl_ultrametric(x),
                                               list(objval = 0)))
    if(.non_ultrametricity(x, max = TRUE) == 0)
        return(.cl_ultrametric_with_meta_added(as.cl_ultrametric(x),
                                               list(objval = 0)))

    ## Handle weights.
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
           SUMT = .l1_fit_ultrametric_by_SUMT(x, weights, control),
           IRIP = .l1_fit_ultrametric_by_IRIP(x, weights, control))
}

### ** .l1_fit_ultrametric_by_SUMT

.l1_fit_ultrametric_by_SUMT <-
function(x, weights = 1, control = list())
{
    ## Try a SUMT with "pseudo-gradients".

    w <- weights / sum(weights)

    ## Control parameters:
    ## gradient,
    gradient <- control$gradient
    if(is.null(gradient))
        gradient <- TRUE
    ## nruns,
    nruns <- control$nruns
    ## start.
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

    ## For the time being, use a simple minimizer.

    n <- attr(x, "Size")
    labels <- attr(x, "Labels")

    L <- function(d) sum(w * abs(d - x))
    P <- .make_penalty_function_ultrametric(n)
    if(gradient) {
        grad_L <- function(d) w * sign(d - x)
        grad_P <- .make_penalty_gradient_ultrametric(n)
    } else
        grad_L <- grad_P <- NULL

    if(is.null(start)) {
        ## Initialize by "random shaking".  Use sd() for simplicity.
        start <- replicate(nruns,  
                           x + rnorm(length(x), sd = sd(x) / sqrt(3)),
                           simplify = FALSE)
    }

    ## And now ...
    out <- sumt(start, L, P, grad_L, grad_P,
                method = control$method, eps = control$eps,
                q = control$q, verbose = control$verbose,
                control = as.list(control$control))

    d <- .ultrametrify(out$x)
    meta <- list(objval = L(d))

    .cl_ultrametric_from_veclh(d, n, labels, meta)
}

### ** .l1_fit_ultrametric_by_IRIP

.l1_fit_ultrametric_by_IRIP <-
function(x, weights = 1, control = list())
{
    ## An attempt of implementing "Iteratively Reweighted Iterative
    ## Projection" as described in Smith (2000, 2001), Journal of
    ## Classification.  Note that this suggests using the Iterative
    ## Projection of Hubert and Arabie (1995), which we cannot as we
    ## have not (yet?) implemented this for the weighted case.  Hence,
    ## we use our SUMT least squares ultrametric fitter instead.
    ##
    ## However, we never got this to converge properly ...
    
    w <- weights / sum(weights)
    
    ## Control parameters:
    ## MIN,
    MIN <- control$MIN
    if(is.null(MIN))
        MIN <- 1e-3
    ## (A rather small cut-off which worked best in the cases we tried.)
    ## eps,
    eps <- control$eps
    if(is.null(eps))
        eps <- 1e-6
    ## maxiter,
    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 100L
    ## reltol,
    reltol <- control$reltol
    if(is.null(reltol))
        reltol <- 1e-6
    ## start,
    start <- control$start
    ## verbose.
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")

    n <- attr(x, "Size")
    labels <- attr(x, "Labels")

    L <- function(d) sum(w * abs(x - d))

    ## Initialize by "random shaking" as for the L2 SUMT, but perhaps we
    ## should not do this?  [Or do it differently?]
    u <- if(is.null(start))
        x + rnorm(length(x), sd = sd(x) / 3)
    else
        start
    ## (No multiple runs for the time being.)
    L_new <- L(u)

    iter <- 1L
    while(iter <= maxiter) {
        if(verbose)
            message(gettextf("Outer iteration: %d", iter))
        L_old <- L_new
        u_old <- u
        weights <- w / pmax(abs(u - x), MIN)
        u <- .ls_fit_ultrametric_by_SUMT(x,
                                         weights = weights,
                                         control = as.list(control$control))
        ## Use some control arguments lateron ...
        L_new <- L(u)
        delta_L <- L_old - L_new
        delta_u <- max(abs(u_old - u))
        if(verbose)
            message(gettextf("Change: u: %g L: %g", delta_u, delta_L))
        if((delta_u < eps) ||
           ((delta_L >= 0) &&
            (delta_L <= reltol * (abs(L_old) + reltol))))
            break
        iter <- iter + 1L
    }

    d <- .ultrametrify(u)
    meta <- list(objval = L(d), status = as.integer(iter == maxiter))

    .cl_ultrametric_from_veclh(d, n, labels, meta)
}

## * ls_fit_sum_of_ultrametrics

ls_fit_sum_of_ultrametrics <-
function(x, nterms = 1, weights = 1, control = list())
{
    if(!inherits(x, "dist"))
        x <- as.dist(x)

    ## We could catch some special cases right away: if x already is an
    ## ultrametric then the fit would be a list with x and nterms - 1
    ## zero ultrametrics ...
    
    ## Control parameters:
    ## eps,
    eps <- control$eps
    if(is.null(eps))
        eps <- 1e-6
    ## maxiter,
    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 100L
    ## method,
    method <- control$method
    if(is.null(method))
        method <- "SUMT"
    ## reltol,
    reltol <- control$reltol
    if(is.null(reltol))
        reltol <- 1e-6
    ## verbose.
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")
    ## Do this at last.
    control <- as.list(control$control)
    ## And be nice ...
    if(identical(method, "SUMT") && is.null(control$nruns))
        control$nruns <- 10L

    L <- function(u)
        sum((x - rowSums(matrix(unlist(u), ncol = nterms))) ^ 2)
    
    ## Init.
    u <- rep.int(list(as.cl_ultrametric(0 * x)), nterms)
    L_new <- L(u)

    ## Loop.
    iter <- 1L
    while(iter <= maxiter) {
        if(verbose)
            message(gettextf("Iteration: %d", iter))
        L_old <- L_new
        delta_u <- 0
        for(i in seq_len(nterms)) {
            if(verbose)
                message(gettextf("Term: %d", i))
            u_old <- u[[i]]
            ## Compute residual r = x - \sum_{j: j \ne i} u(j)
            r <- x - rowSums(matrix(unlist(u[-i]), ncol = nterms - 1L))
            ## Fit residual.
            u[[i]] <- ls_fit_ultrametric(r, method, weights, control)
            ## Accumulate change.
            change <- max(abs(u[[i]] - u_old))
            if(verbose)
                message(gettextf("Change: %g", change))
            delta_u <- max(delta_u, change)
        }
        L_new <- L(u)
        delta_L <- L_old - L_new
        if(verbose)
            message(gettextf("Overall change: u: %g L: %g\n",
                             delta_u, delta_L))
        if((delta_u < eps) ||
           ((delta_L >= 0) &&
            (delta_L <= reltol * (abs(L_old) + reltol))))
            break
        iter <- iter + 1L
    }
    
    .structure(u, objval = L_new, status = as.integer(iter == maxiter))
}

### * as.dist.hclust

## Using hclust() with methods 'median' or 'centroid' typically gives
## reversals and hence not valid hierarchies, i.e., distances which do
## not satisfy the ultrametricity conditions.  The distances can be
## obtained via cophenetic(), but ls_fit_ultrametric() prefers using
## as.dist() [as arguably more appropriate] which in turn can be made to
## "work" by providing as.matrix() methods [bypassing the need to handle
## the extra arguments 'diag' and 'upper' for as.dist()].

as.matrix.hclust <-
function(x, ...)
    as.matrix(cophenetic(x))

### * .non_ultrametricity

.non_ultrametricity <-
function(x, max = FALSE)
{
    if(!is.matrix(x))
        x <- .symmetric_matrix_from_veclh(x)
    .C(C_deviation_from_ultrametricity,
       as.double(x),
       as.integer(nrow(x)),
       fn = double(1L),
       as.logical(max))$fn
}

### * .cl_ultrametric_from_classes

.cl_ultrametric_from_classes <-
function(x)
{
    ## Compute an ultrametric from a hierarchy of classes (i.e., an
    ## n-tree).

    labels <- attr(x, "labels")

    ## Ensure we have no duplicates.
    x <- x[!duplicated(x)]
    
    ## .get_classes_in_hierarchy() orders according to cardinality, but
    ## a consensus method may forget to ...
    x[] <- x[order(sapply(x, length))]

    ## Get the objects (unique codes in the classes).
    objects <- sort(unique(unlist(x)))
    ## (Could also look at the classes of length 1.)

    ## Recursively compute the heights of the classes.
    heights <- double(length = length(x))
    for(i in which(sapply(x, length) > 1)) {
        ## Find the relevant classes.
        j <- sapply(x[seq_len(i - 1L)],
                    function(s) all(s %in% x[[i]]))
        heights[i] <- max(heights[j]) + 1
    }

    ## Next, create an incidence matrix (objects by classes).
    incidences <- sapply(x, function(s) objects %in% s)

    ## Now that we have the heights and incidences, we can compute
    ## distances, using the idea that
    ##     distance(i, j) = min(height(A): A contains i and j)
    n <- length(objects)
    d <- matrix(0, n, n)
    for(i in objects)
        d[i, ] <- heights[apply((rep(incidences[i, ], each = n)
                                 & incidences),
                                1L, which.max)]
    dimnames(d) <- rep.int(list(labels), 2L)

    as.cl_ultrametric(d)
}

### * .cl_ultrametric_with_meta_added

.cl_ultrametric_with_meta_added <-
function(x, meta = NULL)
{
    ## An alternative to adding a 'meta' argument to cl_ultrametric().
    attr(x, "meta") <- meta
    x
}

### .ultrametrify

.ultrametrify <-
function(x)
{
    ## Ensure ultrametricity.
    ## In earlier versions, function
    ## .cl_ultrametric_from_ultrametric_approximation() tried rounding
    ## to non-ultrametric significance, using
    ##   round(x, floor(abs(log10(.non_ultrametricity(x, max = TRUE)))))
    ## which is nice but does not guarantee ultrametricity (and may
    ## result in poorer approximations than what we use now).
    ## Hence, let us use single linkage hierarchical clustering which
    ## gives the best dominated ultrametric approximation.
    cophenetic(hclust(.dist_from_vector(x), "single"))
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
