### * cl_dissimilarity

cl_dissimilarity <-
function(x, y = NULL, method = "euclidean", ...)
{
    x <- as.cl_ensemble(x)
    is_partition_ensemble <-
        (inherits(x, "cl_partition_ensemble")
         || all(sapply(x, .has_object_memberships)))

    ## Be nice.
    if(is.character(y) || is.function(y)) {
        method <- y
        y <- NULL
    }

    if(is.function(method))
        method_name <- "user-defined method"
    else {
        if(!inherits(method, "cl_dissimilarity_method")) {
            ## Get the method definition and description from the
            ## registry. 
            type <- ifelse(is_partition_ensemble,
                           "partition", "hierarchy")
            method <- get_cl_dissimilarity_method(method, type)
        }
        method_name <- method$description
        method <- method$definition
    }
              
    if(!is.null(y)) {
        y <- as.cl_ensemble(y)
        is_partition_ensemble_y <-
            (inherits(y, "cl_partition_ensemble")
             || all(sapply(x, .has_object_memberships)))
        if(!identical(is_partition_ensemble, is_partition_ensemble_y))
            stop("Cannot mix partitions and hierarchies.")
        if(n_of_objects(x) != n_of_objects(y))
            stop("All clusterings must have the same number of objects.")
        ## Build a cross-proximity object of cross-dissimilarities.
        d <- matrix(0, length(x), length(y))
        for(j in seq_along(y))
            d[, j] <- sapply(x, method, y[[j]], ...)
        dimnames(d) <- list(names(x), names(y))
        return(cl_cross_proximity(d, method_name,
                                  class = "cl_cross_dissimilarity"))
    }
    
    ## Otherwise, build a proximity object of dissimilarities.
    n <- length(x)
    d <- vector("list", length = n - 1L)
    ind <- seq_len(n)
    while(length(ind) > 1L) {
        j <- ind[1L]
        ind <- ind[-1L]
        d[[j]] <- sapply(x[ind], method, x[[j]], ...)
    }
    
    cl_proximity(unlist(d),
                 method_name,
                 labels = names(x),
                 size = n,
                 class = c("cl_dissimilarity", "cl_proximity", "dist"))
}

### ** .cl_dissimilarity_partition_euclidean

.cl_dissimilarity_partition_euclidean <-
function(x, y)
{
    k <- max(n_of_classes(x), n_of_classes(y))
    M_x <- cl_membership(x, k)
    M_y <- cl_membership(y, k)
    ## Match classes from conforming memberships.
    ind <- solve_LSAP(crossprod(M_x, M_y), maximum = TRUE)
    sqrt(sum((M_x - M_y[, ind]) ^ 2))
}

### ### ** .cl_dissimilarity_partition_manhattan

.cl_dissimilarity_partition_manhattan <-
function(x, y)
{
    k <- max(n_of_classes(x), n_of_classes(y))
    M_x <- cl_membership(x, k)
    M_y <- cl_membership(y, k)
    C <- .cxdist(M_x, M_y, "manhattan")
    ind <- solve_LSAP(C)
    sum(C[cbind(seq_along(ind), ind)])
}    

### ** .cl_dissimilarity_partition_comemberships

.cl_dissimilarity_partition_comemberships <-
function(x, y)
{
    ## We used to have the straightforward
    ##   C_x <- tcrossprod(cl_membership(x)) # M_x M_x'
    ##   C_y <- tcrossprod(cl_membership(y)) # M_y M_y'
    ##   sum((C_x - C_y) ^ 2) / n_of_objects(x) ^ 2
    ## But note that
    ##   \| AA' - BB' \|^2
    ##      = tr((AA' - BB')'(AA' - BB')
    ##      = tr(A'A A'A) - 2 tr(A'B B'A) + tr(B'B B'B)
    ##      = \| A'A \|^2 - 2 \| A'B \|^2 + \| B'B \|^2
    ## which can be computed much more efficiently as all involved cross
    ## product matrices are "small" ...
    k <- max(n_of_classes(x), n_of_classes(y))
    M_x <- cl_membership(x, k)
    M_y <- cl_membership(y, k)
    sqrt(sum(crossprod(M_x) ^ 2)
         - 2 * sum(crossprod(M_x, M_y) ^ 2)
         + sum(crossprod(M_y) ^ 2))
}

### ** .cl_dissimilarity_partition_symdiff

.cl_dissimilarity_partition_symdiff <-
function(x, y)
{
    ## Cardinality of the symmetric difference of the partitions
    ## regarded as binary equivalence relations, i.e., the number of
    ## discordant pairs.

    ## Handle soft partitions using the corresponding hard ones.
    ## (At least, for the time being.)

    n <- n_of_objects(x)
    .cl_dissimilarity_partition_Rand(x, y) * choose(n, 2)
}
    
### ** .cl_dissimilarity_partition_Rand

.cl_dissimilarity_partition_Rand <-
function(x, y)    
{
    ## Handle soft partitions using the corresponding hard ones.
    ## (At least, for the time being.)

    1 - .cl_agreement_partition_Rand(x, y)
}

### ** .cl_dissimilarity_partition_GV1

.cl_dissimilarity_partition_GV1 <-
function(x, y)
{    
    k_x <- n_of_classes(x)
    k_y <- n_of_classes(y)
    M_x <- cl_membership(x, k_x)
    M_y <- cl_membership(y, k_y)
    C <- outer(colSums(M_x ^ 2), colSums(M_y ^ 2), "+") -
        2 * crossprod(M_x, M_y)
    if(k_x < k_y)
        C <- rbind(C, matrix(0, nrow = k_y - k_x, ncol = k_y))
    else if(k_x > k_y)
        C <- cbind(C, matrix(0, nrow = k_x, ncol = k_x - k_y))
    ind <- solve_LSAP(C)
    sqrt(sum(C[cbind(seq_along(ind), ind)]))
    ## (Note that this sum really only includes matched non-dummy
    ## classes.)
}

### ** .cl_dissimilarity_partition_BA_A

.cl_dissimilarity_partition_BA_A <-
function(x, y)
{
    .cl_dissimilarity_partition_manhattan(as.cl_hard_partition(x),
                                          as.cl_hard_partition(y)) / 2
    ## Could to this more efficiently, of course ...
}

### ** .cl_dissimilarity_partition_BA_C

.cl_dissimilarity_partition_BA_C <-
function(x, y)
{
    n_of_classes(x) + n_of_classes(y) - 2 * n_of_classes(cl_join(x, y))
}

### ** .cl_dissimilarity_partition_BA_D

.cl_dissimilarity_partition_BA_D <-
    .cl_dissimilarity_partition_Rand

### ** .cl_dissimilarity_partition_BA_E

.cl_dissimilarity_partition_BA_E <-
function(x, y)
{
    z <- table(cl_class_ids(x), cl_class_ids(y))
    z <- z / sum(z)

    ## Average mutual information between the partitions.
    y <- outer(rowSums(z), colSums(z))
    i <- which((z > 0) & (y > 0))
    I <- sum(z[i] * log(z[i] / y[i]))

    ## Entropy of meet(x, y).
    i <- which(z > 0)
    H <- - sum(z[i] * log(z[i]))

    1 - I / H
}

### ** .cl_dissimilarity_partition_VI

.cl_dissimilarity_partition_VI <-
function(x, y, weights = 1)
{
    ## Variation of information for general "soft clusterings", cf
    ## Section 5.2. in Meila (2002). 
    weights <- rep(weights, length.out = n_of_objects(x))
    weights <- weights / sum(weights)
    M_x <- cl_membership(x)
    ## Weighted marginal distribution of x:
    m_x <- colSums(weights * M_x)
    M_y <- cl_membership(y)
    ## Weighted marginal distribution of y:
    m_y <- colSums(weights * M_y)
    gamma <- crossprod(weights * M_x, M_y)
    delta <- outer(m_x, m_y)
    ## Entropy of x:
    H_x <- - sum(m_x * log(ifelse(m_x > 0, m_x, 1)))
    ## Entropy of y:
    H_y <- - sum(m_y * log(ifelse(m_y > 0, m_y, 1)))
    ## VI is H_x + H_y minus twice the (weighted) joint information.
    i <- which((gamma > 0) & (delta > 0))
    H_x + H_y - 2 * sum(gamma[i] * log(gamma[i] / delta[i]))
}

### ** .cl_dissimilarity_partition_Mallows

.cl_dissimilarity_partition_Mallows <-
function(x, y, p = 1, alpha = NULL, beta = NULL)
{
    ## Currently, no "real" primal-dual solver for minimum cost flow
    ## problems, and lpSolve::lp.transport() seems to work only for
    ## integer bounds.  Hence, rather than using
    ##
    ##   C <- .cxdist(cl_membership(x), cl_membership(y),
    ##                "minkowski", p) ^ p
    ##   n_x <- nrow(C)
    ##   n_y <- ncol(C)
    ##   if(is.null(alpha))
    ##       alpha <- rep.int(1 / n_x, n_x)
    ##   else {
    ##       alpha <- rep(alpha, length.out = n_x)
    ##       alpha <- alpha / sum(alpha)
    ##   }
    ##
    ## etc right away, ensure a square cost matrix so that we can have
    ## integer bounds for at least the default case.

    k <- max(n_of_classes(x), n_of_classes(y))
    M_x <- cl_membership(x, k)
    M_y <- cl_membership(y, k)
    C <- .cxdist(M_x, M_y, "minkowski", p) ^ p
    if(is.null(alpha)) alpha <- rep.int(1, k)
    if(is.null(beta)) beta <- rep.int(1, k)
    lpSolve::lp.transport(C, "min",
                          rep("==", k), alpha,
                          rep("==", k), beta,
                          integers = NULL)$objval ^ (1 / p)
}

### ** .cl_dissimilarity_partition_CSSD

.cl_dissimilarity_partition_CSSD <-
function(x, y, L = NULL, alpha = NULL, beta = NULL, ...)
{
    ## Cluster Similarity Sensitive Distance.
    ## Reference:  D. Zhou, J. Li and H. Zha (2005),
    ##  A new Mallows distance based metric for comparing clusterings.
    ## See .cl_dissimilarity_partition_Mallows() re solving cost flow
    ## problems.

    ## Dissimilarity is defined by minimizing
    ##   \sum_{k,l} (1 - 2 w_{kl} / (alpha_k + beta_l)) L_{kl}
    ## where
    ##   L_{kl} = \sum_i m_{x;ik} m_{y;il} distance(p_{x;k}, p_{y;l})
    ## with m and p the memberships and prototypes, respectively.
    ## If we get matrices of prototypes, use .rxdist; otherwise, the
    ## user needs to specify an L function or matrix.

    k_x <- n_of_classes(x)
    k_y <- n_of_classes(y)
    M_x <- cl_membership(x, k_x)
    M_y <- cl_membership(y, k_y)
    if(!is.matrix(L)) {
        p_x <- cl_prototypes(x)
        p_y <- cl_prototypes(y)
        if(is.matrix(p_x) && is.matrix(p_y) && is.null(L))
            L <- .rxdist(p_x, p_y, ...)
        else if(is.function(L))
            L <- L(p_x, p_y)
        else
            stop("Cannot compute prototype distances.")
    }
    C <- crossprod(M_x, M_y) * L
    if(is.null(alpha)) alpha <- rep.int(1, k_x)
    if(is.null(beta)) beta <- rep.int(1, k_y)
    sum(C) - 2 * lpSolve::lp.transport(C / outer(alpha, beta, "+"),
                                       "max",
                                       rep("==", k_x), alpha,
                                       rep("==", k_y), beta,
                                       integers = NULL)$objval
}

### ** .cl_dissimilarity_hierarchy_euclidean

.cl_dissimilarity_hierarchy_euclidean <-
function(x, y, weights = 1)
{
    if(!.has_object_dissimilarities(x) ||
       !.has_object_dissimilarities(y))
        return(NA)
    u <- cl_object_dissimilarities(x)
    v <- cl_object_dissimilarities(y)
    sqrt(sum(weights * (u - v) ^ 2))
}

### ** .cl_dissimilarity_hierarchy_manhattan

.cl_dissimilarity_hierarchy_manhattan <-
function(x, y, weights = 1)
{
    if(!.has_object_dissimilarities(x) ||
       !.has_object_dissimilarities(y))
        return(NA)
    u <- cl_object_dissimilarities(x)
    v <- cl_object_dissimilarities(y)
    sum(weights * abs(u - v))
}

### ** .cl_dissimilarity_hierarchy_cophenetic

.cl_dissimilarity_hierarchy_cophenetic <-
function(x, y)
{
    if(!.has_object_dissimilarities(x) ||
       !.has_object_dissimilarities(y))
        return(NA)
    u <- cl_object_dissimilarities(x)
    v <- cl_object_dissimilarities(y)
    1 - cor(u, v) ^ 2
}

### ** .cl_dissimilarity_hierarchy_gamma

.cl_dissimilarity_hierarchy_gamma <-
function(x, y)
{
    ## <NOTE>
    ## This is a dissimilarity measure that works for arbitrary
    ## dissimilarities, see e.g. Bock.
    ## (And the current implementation finally respects this ...)
    ## </NOTE>
    if(!.has_object_dissimilarities(x) ||
       !.has_object_dissimilarities(y))
        return(NA)
    u <- cl_object_dissimilarities(x)
    v <- cl_object_dissimilarities(y)
    n <- length(u)
    .C(C_clue_dissimilarity_count_inversions,
       as.double(u),
       as.double(v),
       as.integer(n),
       count = double(1L)) $ count / choose(n, 2)
}

### ** .cl_dissimilarity_hierarchy_symdiff

.cl_dissimilarity_hierarchy_symdiff <-
function(x, y)
{
    ## Cardinality of the symmetric difference of the n-trees when
    ## regarded as sets of subsets (classes) of the set of objects.

    x <- cl_classes(x)
    y <- cl_classes(y)
    sum(is.na(match(x, y))) + sum(is.na(match(y, x)))
}

### ** .cl_dissimilarity_hierarchy_Chebyshev

.cl_dissimilarity_hierarchy_Chebyshev <-
function(x, y)
{
    if(!.has_object_dissimilarities(x) ||
       !.has_object_dissimilarities(y))
        return(NA)
    u <- cl_object_dissimilarities(x)
    v <- cl_object_dissimilarities(y)
    max(abs(u - v))
}

### ** .cl_dissimilarity_hierarchy_Lyapunov

.cl_dissimilarity_hierarchy_Lyapunov <-
function(x, y)
{
    if(!.has_object_dissimilarities(x) ||
       !.has_object_dissimilarities(y))
        return(NA)
    q <- cl_object_dissimilarities(x) / cl_object_dissimilarities(y)
    if(is.matrix(q)) q <- q[lower.tri(q)]
    log(max(q) / min(q))
}

### ** .cl_dissimilarity_hierarchy_BO

.cl_dissimilarity_hierarchy_BO <-
function(x, y, delta, ...)
{
    ## Compute Boorman-Olivier (1973) dendrogram ("valued tree")
    ## dissimilarities of the form
    ##
    ##    m_\delta(T_1, T_2)
    ##      = \int_0^\infty \delta(P_1(\alpha), P_2(\alpha)) d\alpha
    ##
    ## where the trees (dendrograms) are defined as right-continuous
    ## maps from [0, \Infty) to the partition lattice.

    ## We can compute this as follows.  Take the ultrametrics and use
    ## as.hclust() to detemine the heights \alpha_1(k) and \alpha_2(l)
    ## of the splits.  Let \alpha_i be the sequence obtained by
    ## combining these two.  Then
    ##
    ##   m_\delta
    ##     = \sum_{i=0}^{L-1} (\alpha_{i+1} - \alpha_i)
    ##                        \delta(P_1(\alpha_i), P_2(\alpha_i))
    ##
    ## We use cutree() for computing the latter partitions.  As we
    ## already have the hclust representations, we should be able to do
    ## things more efficiently ...

    if(inherits(x, "hclust"))
        t_x <- x
    else if(inherits(x, "cl_ultrametric"))
        t_x <- as.hclust(x)
    else if(is.cl_dendrogram(x))
        t_x <- as.hclust(cl_ultrametric(x))
    else
        return(NA)
    if(inherits(y, "hclust"))
        t_y <- y
    else if(inherits(y, "cl_ultrametric"))
        t_y <- as.hclust(y)
    else if(is.cl_dendrogram(y))
        t_y <- as.hclust(cl_ultrametric(y))
    else
        return(NA)

    if(is.unsorted(t_x$height) || is.unsorted(t_y$height))
        return(NA)
    
    alpha <- sort(unique(c(t_x$height, t_y$height)))
    cuts_x <- cutree(t_x, h = alpha)
    cuts_y <- cutree(t_y, h = alpha)
    deltas <- mapply(cl_dissimilarity,
                     lapply(split(cuts_x, col(cuts_x)),
                            as.cl_partition),
                     lapply(split(cuts_y, col(cuts_y)),
                            as.cl_partition),
                     MoreArgs = list(delta, ...))
    sum(diff(alpha) * deltas[-length(deltas)])
}

### ** .cl_dissimilarity_hierarchy_spectral

.cl_dissimilarity_hierarchy_spectral <-
function(x, y)
{
    if(!.has_object_dissimilarities(x) ||
       !.has_object_dissimilarities(y))
        return(NA)
    u <- cl_object_dissimilarities(x)
    v <- cl_object_dissimilarities(y)
    svd(as.matrix(u - v))$d[1L]
}


### * as.dist.cl_dissimilarity

as.dist.cl_dissimilarity <-
function(m, diag = FALSE, upper = FALSE)
{
    y <- c(m)
    ## Fill non-inherited attributes with default values.
    attributes(y) <- c(attributes(m)[c("Size", "Labels")],
                       Diag = diag, Upper = upper, call = match.call())
    ## (Note that as.dist.default() does not automatically add
    ## 'method'.)
    class(y) <- "dist"
    y
}

### * [.cl_dissimilarity

"[.cl_dissimilarity" <-
function(x, i, j)
{
    y <- NextMethod("[")
    if(!inherits(y, "cl_dissimilarity")) {
        description <- attr(x, "description")
        return(cl_cross_proximity(y, description = description,
                                  class = "cl_cross_dissimilarity"))
    }
    y
}

### .cxdist

.cxdist <-
function(A, B, method = c("euclidean", "manhattan", "minkowski"), ...)
{
    ## Return the column cross distance matrix of A and B.
    ## I.e., the matrix C = [c_{j,k}] with
    ##   c_{j,k} = distance(A[, j], B[, k])
    ## Currently, only Manhattan (L1) distances are provided.
    ## Extensions to Minkowski or even more distances (a la dist())
    ## could be added eventually.

    ## <NOTE>
    ## Possible implementations include
    ##
    ## foo_a <- function(A, B)
    ##   apply(B, 2, function(u) colSums(abs(A - u)))
    ## foo_d <- function(A, B) {
    ##   out <- as.matrix(dist(rbind(t(A), t(B)), "manhattan"))
    ##   dimnames(out) <- NULL
    ##   nc_B <- NCOL(B)
    ##   out[seq(from = NCOL(A) + 1, length.out = nc_B), seq_len(nc_B)]
    ## }
    ## foo_f <- function(A, B) {
    ##   out <- matrix(0, NCOL(A), NCOL(B))
    ##   for(j in seq_len(NCOL(A)))
    ##     for(k in seq_len(NCOL(B)))
    ##       out[j, k] = sum(abs(A[, j] - B[, k]))
    ##   out
    ## }
    ##
    ## The one actually used seems to be the best performer, with the
    ## "for" version a close second (note that "typically", A and B have
    ## much fewer columns than rows).
    ## only few columns

    method <- match.arg(method)
    
    ## Workhorse.
    FOO <- switch(method,
                  "euclidean" =
                  function(M) sqrt(colSums(M ^ 2)),
                  "manhattan" = 
                  function(M) colSums(abs(M)),
                  "minkowski" = {
                      ## Power needs to be given.
                      p <- list(...)[[1L]]
                      function(M)
                          (colSums(abs(M) ^ p)) ^ (1 / p)
                  })
    
    
    out <- matrix(0, NCOL(A), NCOL(B))
    for(k in seq_len(NCOL(B)))
        out[, k] <- FOO(A - B[, k])
    out
}

### .rxdist

.rxdist <-
function(A, B, method = c("euclidean", "manhattan", "minkowski"), ...)
{
    ## Return the row cross distance matrix of A and B.
    ## I.e., the matrix C = [c_{j,k}] with
    ##   c_{j,k} = distance(A[j, ], B[k, ])
    
    ## <NOTE>
    ## Could also do something like
    ##   ind <- seq_len(NROW(B))
    ##   as.matrix(dist(rbind(B, A)))[-ind, ind]
    ## but that is *very* inefficient for the "usual" data by prototype
    ## case (where NROW(B) << NROW(A)).
    ## </NOTE>
    
    ## No fancy pmatching for methods for the time being.
    method <- match.arg(method)
    
    ## Workhorse: Full A, single row of b.
    FOO <- switch(method,
                  "euclidean" =
                  function(A, b) sqrt(rowSums(sweep(A, 2, b) ^ 2)),
                  "manhattan" = 
                  function(A, b) rowSums(abs(sweep(A, 2, b))),
                  "minkowski" = {
                      ## Power needs to be given.
                      p <- list(...)[[1L]]
                      function(A, b)
                          (rowSums(abs(sweep(A, 2, b)) ^ p)) ^ (1 / p)
                  })
                      
    out <- matrix(0, NROW(A), NROW(B))
    for(k in seq_len(NROW(B)))
        out[, k] <- FOO(A, B[k, ])
    out
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
