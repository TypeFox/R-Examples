## Spherical (possibly sparse) k-means.
## Authors: Kurt Hornik, Ingo Feinerer, Martin Kober.

## Partition given vectors x_b by minimizing the criterion
##   \sum w_b u_{bj}^m d(x_b, p_j)
## where the w_b are case weights, u_{bj} is the membership of x_b to
## class j, p_j is the prototype of class j (and thus minimizes
##   \sum_b w_b u_{bj}^m d(x_b, p)
## over p), and d is *cosine dissimilarity*
##   d(x, p) = 1 - cos(x, p),
## where
##   cos(x, p) = <x, p> / (||x|| * ||p||)
## (so that d(x, p) is half the squared euclidean distance between the
## normalized x and the normalized p).
##
## (E.g., http://www.ams.jhu.edu/~priebe/.FILES/applepres.pdf or
## http://reference.wolfram.com/mathematica/ref/CosineDistance.html for
## similar use of cosine dissimilarity/distance.)
##
## The literature only considers the unweighted hard spherical k-means
## problem (all w_b = 1 and m = 1), which we refer to as the "standard
## spherical k-means problem".  This has the criterion function
##   \sum_j \sum_{b in class j} d(x_b, p_j)
##     = B - \sum_j \sum_{b in class j} cos(x_b, p_j)
## with B the number of objects, so that minimizing the dissimilarity
## criterion is equivalent to maximizing
##   \sum_j \sum_{b in class j} cos(x_b, p_j))
## which is criterion I2 in CLUTO.
##
## All built-in standard spherical k-means methods internally use I2 but
## report the value of the corresponding dissimilarity criterion.

### * Infrastructure

## Prototypes are always represented by a dense matrix.
## Objects can be represented by dense matrices or sparse matrices
## provided that the following matrix computations work for these:
## * Subscripting rows
## * Scaling rows via multiplication (from the left) or division (from
##   the right) by a numeric vector
## * Taking (element-wise) squares
## * Coercion to dense matrix via as.matrix()
## * rowSums and colSums
## * tcrossprod with a dense matrix.

## Unfortunately, whereas the first four can easily be accomplished via
## (S3) subscript, Ops group and as.matrix() methods, rowSums, colSums
## and tcrossprod are not generic, and S4 dispatch will not occur on the
## base function.  Hence, we provide our own S3 generics and methods for
## supported classes.

## For the CLUTO interface, an as.simple_triplet_matrix() method must be
## available.

g_tcrossprod <-
function(x, y = NULL)
    UseMethod("g_tcrossprod")
g_tcrossprod.default <-
function(x, y = NULL)
    base::tcrossprod(x, y)
g_tcrossprod.simple_triplet_matrix <-
function(x, y = NULL)
    tcrossprod_simple_triplet_matrix(x, y)
g_tcrossprod.dgCMatrix <-
function(x, y = NULL)
    Matrix::tcrossprod(x, y)
g_tcrossprod.dgTMatrix <-
function(x, y = NULL)
    Matrix::tcrossprod(x, y)

## (Using
##    g_tcrossprod.default <- base::tcrossprod
## has R CMD check NOTE .Internal calls, using
##    g_tcrossprod.dgTMatrix <- Matrix::tcrossprod
## loads Matrix dependencies.)

g_row_sums <-
function(x, na.rm = FALSE, dims = 1, ...)
    UseMethod("g_row_sums")
## <NOTE>
## These also also provided by package slam now, so one could simply do
##   g_row_sums.default <-
##   function(x, na.rm = FALSE, dims = 1, ...)
##     slam::row_sums(x, na.rm = na.rm, dims = dims, ...)
g_row_sums.default <-
function(x, na.rm = FALSE, dims = 1, ...)
    base::rowSums(x, na.rm = na.rm, dims = dims, ...)
g_row_sums.simple_triplet_matrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    slam::row_sums(x, na.rm = na.rm, dims = dims, ...)
g_row_sums.dgCMatrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    Matrix::rowSums(x, na.rm = na.rm, dims = dims, ...)
g_row_sums.dgTMatrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    Matrix::rowSums(x, na.rm = na.rm, dims = dims, ...)
## </NOTE>

g_col_sums <-
function(x, na.rm = FALSE, dims = 1, ...)
    UseMethod("g_col_sums")
## <NOTE>
## These also also provided by package slam now, so one could simply do
##   g_col_sums.default <-
##   function(x, na.rm = FALSE, dims = 1, ...)
##     slam::col_sums(x, na.rm = na.rm, dims = dims, ...)
g_col_sums.default <-
function(x, na.rm = FALSE, dims = 1, ...)
    base::colSums(x, na.rm = na.rm, dims = dims, ...)
g_col_sums.simple_triplet_matrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    slam::col_sums(x, na.rm = na.rm, dims = dims, ...)
g_col_sums.dgCMatrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    Matrix::colSums(x, na.rm = na.rm, dims = dims, ...)
g_col_sums.dgTMatrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    Matrix::colSums(x, na.rm = na.rm, dims = dims, ...)
## </NOTE>

g_col_sums_by_group <-
function(x, g, ...)
    UseMethod("g_col_sums_by_group")
g_col_sums_by_group.simple_triplet_matrix <-
function(x, g, ...)
    as.matrix(rollup(x, MARGIN = 1L, INDEX = factor(g), FUN = sum,
                     na.rm = FALSE))
g_col_sums_by_group.default <-
function(x, g, ...)
{
    g <- factor(g)
    v <- levels(g)
    out <- matrix(0, length(v), ncol(x), dimnames = list(v, colnames(x)))
    for(i in seq_along(v))
        out[i, ] <- g_col_sums_with_logical_index(x, g == v[i])
    out
}

g_crossprod <-
function(x, y = NULL)
    UseMethod("g_crossprod", y)
g_crossprod.default <-
function(x, y = NULL)
    base::crossprod(x, y)
g_crossprod.simple_triplet_matrix <-
function(x, y = NULL)
    slam::crossprod_simple_triplet_matrix(x, y)
g_crossprod.dgCMatrix <-
function(x, y = NULL)
    Matrix::crossprod(x, y)
g_crossprod.dgTMatrix <-
function(x, y = NULL)
    Matrix::crossprod(x, y)

### * skmeans_xdist

skmeans_xdist <-
function(x, y = NULL)
{
    x <- row_normalize(x)
    if(!is.null(y))
        y <- row_normalize(y)
    pmax(1 - g_tcrossprod(x, y), 0)
}

### * skmeans_family

## In the skmeans family object, we provide the "correct" cosine
## dissimilarity between objects and prototypes, irrespective of
## possible normalization.  For the methods, we internally always
## normalize objects and prototypes and hence use faster D and C
## functions.

skmeans_family <-
    pclust_family(D =
                  function(x, prototypes) {
                      skmeans_xdist(x, prototypes)
                  },
                  C =
                  function(x, weights, control) {
                      ## Computes weighted averages of the
                      ## normalized data.
                      x <- row_normalize(x)
                      weights <- weights / sum(weights)
                      out <- g_col_sums_with_weights(x, weights)
                      ## <NOTE>
                      ## Should prototypes be normalized?
                      ## They are not in the basic references.
                      ## If normalization is desired, uncomment
                      ##    out <- out / sqrt(sum(out ^ 2))
                      ## </NOTE>
                      out
                  },
                  init =
                  function(x, k) {
                      ## <NOTE>
                      ## Should perhaps ensure that this returns
                      ## unique prototypes.
                      as.matrix(x[sample.int(nrow(x), k), ,
                                  drop = FALSE])
                      ## </NOTE>
                  },
                  description = "spherical k-means")

### * skmeans

skmeans <-
function(x, k, method = NULL, m = 1, weights = 1, control = list())
{
    mc <- match.call()
    
    if(!all(row_norms(x) > 0))
        stop("Zero rows are not allowed.")

    ## Methods must at least have formals x, k and control.
    ## Try to ensure that formals m and weights are only used if
    ## different from the default ("simple") case.

    args <- list(x = x, k = k, control = control)
    if(!missing(m) && !identical(m, 1))
        args <- c(args, list(m = m))
    if(!missing(weights) && !all(weights == 1))
        args <- c(args,
                  list(weights = rep(weights, length.out = nrow(x))))

    skmeans_methods <-
        c(genetic = "genetic",
          pclust = "pclust",
          CLUTO = "CLUTO",
          gmeans = "gmeans",
          kmndirs = "kmndirs",
          ## <FIXME>
          ## Remove eventually.
          LIH = "local_improvement_heuristic",
          LIHC = "local_improvement_heuristic_with_chains"
          ## </FIXME>
          )

    if(!is.function(method)) {
        method <- if(is.null(method)) {
            if(is.null(args$m)) "genetic" else "pclust"
        } else if(is.character(method)) {
            pos <- pmatch(tolower(method),
                          tolower(names(skmeans_methods)))
            if(is.na(pos))
                stop(gettextf("Invalid skmeans method '%s'.", method),
                     domain = NA)
            method <- skmeans_methods[pos]
        } else {
            stop("Invalid skmeans method.")
        }
        method <- get(sprintf(".skmeans_%s", method))
    }

    ## Now check whether the method has the required formals.
    na <- names(args)
    nf <- names(formals(method))
    if(any(ind <- is.na(match(na, nf))))
        stop(gettextf("Given skmeans method lacks formals %s",
                      paste(sQuote(na[ind]), collapse = " and ")),
             domain = NA)
    ## Otherwise, ensure that method gets called with "all" arguments so
    ## that it does not also have to provide defaults for these.
    if(("m" %in% nf) && !("m" %in% na))
        args <- c(args, list(m = m))
    if(("weights" %in% nf) && !("weights" %in% na))
        args <- c(args,
                  list(weights = rep(weights, length.out = nrow(x))))

    ## Call the skmeans method.
    y <- do.call(method, args)
    ## Record the original call.
    y$call <- mc

    y
}

### * skmeans methods

print.skmeans <-
function(x, ...)
{
    ids <- x$cluster
    m <- x$m
    tab <- table(ids)
    sizes <- paste(tab, collapse = ", ")
    txt <- if(m == 1) {
        c(gettextf("A hard spherical k-means partition of %d objects into %d classes.", 
                 length(ids), length(tab)),
          gettextf("Class sizes: %s", sizes))
    } else {
        c(gettextf("A soft spherical k-means partition (degree m = %f) of %d objects into %d classes.",
                   m, length(ids), length(tab)),
          gettextf("Class sizes of closest hard partition: %s.",
                   sizes))
    }
    writeLines(strwrap(txt))
    cat("Call: ", paste(deparse(x$call), collapse = "\n"), "\n",
        sep = "")
    invisible(x)
}

cl_validity.skmeans <-
function(x, data = NULL, ...)
{
    ## Dissimilarity accounted for by spherical k-means partitions.

    ## If the original data matrix is not given, try to obtain it from
    ## the orginal skmeans() call.
    if(is.null(data))
        data <- eval(x$call$x, parent.frame())

    data <- row_normalize(data)

    n <- nrow(data)

    v <- if(x$m == 1) {
        ## Hard partitions.
        ids <- factor(x$cluster)
        sizes <- table(ids)
        k <- length(sizes)

        sums <- g_col_sums_by_group(data, ids)
        
        1 - ((1 - sum(sums ^ 2) / sum(sizes ^ 2)) /
             (1 - sum(g_col_sums(data) ^ 2) / n ^ 2))
    } else {
        ## Soft partitions.
        M <- cl_membership(x)

        sums <- g_crossprod(M, data)

        1 - ((1 - sum(sums ^ 2) / sum(colSums(M) ^ 2)) /
             (1 - sum(g_col_sums(data) ^ 2) / n ^ 2))
    }

    out <- list("Dissimilarity accounted for" = v)
    class(out) <- "cl_validity"

    out
}

silhouette.skmeans <-
function(x, data = NULL, ...)
{
    mc <- match.call()

    ## If the original data matrix is not given, try to obtain it from
    ## the orginal skmeans() call.
    if(is.null(data))
        data <- eval(x$call$x, parent.frame())

    data <- row_normalize(data)

    n <- nrow(data)

    ids <- factor(x$cluster)
    sizes <- table(ids)
    k <- length(sizes)

    sums <- g_col_sums_by_group(data, ids)

    ind <- cbind(seq_len(n), ids)
    numers <- denoms <- matrix(rep.int(sizes, rep.int(n, k)), n, k)
    denoms[ind] <- sizes[ids] - 1L

    aves <- (numers - g_tcrossprod(data, sums)) / denoms

    a <- aves[ind]
    aves[ind] <- Inf
    neighbor <- max.col(-aves)
    b <- aves[cbind(seq_len(n), neighbor)]
    s <- (b - a) / pmax(a, b)
    s[is.nan(s)] <- 0

    wds <- cbind(cluster = ids, niehgbor = neighbor, sil_width = s)
    rownames(wds) <- names(ids)
    attr(wds, "Ordered") <- FALSE
    attr(wds, "call") <- mc
    class(wds) <- "silhouette"

    wds
}

### * .skmeans_pclust

.skmeans_pclust <-
function(x, k, m = 1, weights = 1, control = NULL)
{
    x <- row_normalize(x)
    
    ## Handle starting values.
    start <- control$start
    if(is.null(start)) {
        nruns <- control$nruns
        if(is.null(nruns))
            nruns <- 1L
        start <- as.list(rep.int("p", nruns))
    } else if(!is.list(start) || inherits(start, "skmeans"))
        start <- list(start)
    
    control$start <- start
    control$nruns <- NULL

    if(m == 1)
        .skmeans_hard_pclust(x, k, weights, control)
    else
        .skmeans_soft_pclust(x, k, m, weights, control)
}

### * .skmeans_hard_pclust

.skmeans_hard_pclust <-
function(x, k, weights = 1, control = NULL)
{
    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 100L
    
    maxchains <- control$maxchains
    if(is.null(maxchains))
        maxchains <- 0L
    
    reltol <- control$reltol
    if(is.null(reltol))
        reltol <- sqrt(.Machine$double.eps)
    
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")

    nr <- nrow(x)
    nc <- ncol(x)

    ## In the weighted case, we can perform all computations on w_i x_i
    ## (after normalization).  For first variation moves, we need w_i^2.
    if(all(weights == 1)) {
        wsq <- rep.int(1, nr)
        s <- nr
    } else {
        if(!all(weights > 0))
            stop("All weights must be positive.")
        wsq <- weights ^ 2
        s <- sum(weights)
        x <- weights * x
    }
    
    ## Initialize.
    start <- control$start
    nruns <- length(start)
    
    opt_value <- 0
    run <- 1L

    if(verbose && (nruns > 1L))
        message(gettextf("Pclust run: %d", run))

    repeat {
        p <- .skmeans_init_for_normalized_x(x, k, start[[run]], weights)
        old_value <- 0
        iter <- 1L
        while(iter <= maxiter) {
            ## Fixed point iteration.
            similarities <- g_tcrossprod(x, p)
            ## New ids.
            ids <- ids_from_similarities(similarities, k)
            ## New prototypes.
            sums <- .hard_skmeans_C_for_normalized_x(x, ids)
            norms <- row_norms(sums)
            p <- sums / norms
            ## New value.
            new_value <- sum(norms)
            
            if(verbose)
                message(gettextf("Iteration: %d *** value: %g",
                                 iter, s - new_value))

            ## If the change from the block update was big enough,
            ## continue with the next fixed point iteration.
            ## Otherwise, if maxchains is positive, try first variation
            ## steps.
            if(abs(old_value - new_value)
               >= reltol * (abs(old_value) + reltol)) {
                old_value <- new_value                
                iter <- iter + 1L
                next
            }
            else if(maxchains == 0L)
                break

            ## Try first variation steps.
            ## What we need to determine is the optimal chain of first
            ## variation moves of length up to maxchains.  If this chain
            ## results in a change "big enough" (above the tolerance) we
            ## perform the move and resume fixed point iteration;
            ## otherwise, we are done.
            ## We do this as follows:
            ## * Iterate over chains from 1 to maxchains.
            ## * For chains giving a change which is "big enough": if
            ##   this change is the currenly optimal one, record the
            ##   chain number (opt_chains) and update ids, p and norms.
            ## * At the end: if no optimal change was recorded, we are
            ##   done.
            opt_chains <- 0L
            opt_change <- 0
            tol <- reltol * (abs(new_value) + reltol)
            fv_ids <- ids
            fv_sums <- sums
            fv_norms <- norms
            fv_p <- p
            crossprods <- 2 * sweep(similarities, 2L, norms, "*")
            ## Could special-case maxchains == 1 to save a few msecs,
            ## but prefer clarity for speed here.
            change <- 0            
            marked <- integer()
            chains <- 1L
            while(chains <= maxchains) {
                nids <- fv_norms[fv_ids]
                ind <- cbind(seq_len(nr), fv_ids)
                ## The effects of removing an object from its cluster.
                new_norms_rem <-
                    sqrt(nids ^ 2 - crossprods[ind] + wsq[fv_ids])
                Delta_Q_rem <- new_norms_rem - nids
                ## The effects of adding an object to another cluster.
                new_norms_add <-
                    sqrt(outer(wsq, fv_norms ^ 2, "+") + crossprods)
                Delta_Q_add <- sweep(new_norms_add, 2L, fv_norms, "-")
                ## What is the best such move?
                Delta_Q <- Delta_Q_rem + Delta_Q_add
                Delta_Q[ind] <- -1e9
                Delta_Q[marked, ] <- -1e9
                Delta_Q <- as.numeric(Delta_Q)
                pos <- which.max(Delta_Q)
                ## Change object o from its cluster i (fv_ids[o]) to
                ## cluster j.
                change <- Delta_Q[pos] + change
                j <- (pos - 1L) %/% nr + 1L
                o <- (pos - 1L) %% nr + 1L
                i <- fv_ids[o]
                marked <- c(marked, o)
                if(verbose)
                    message(c(gettextf("Iteration: %d [CH %d] *** value: %g change: %g",
                                       iter, chains, s - new_value,
                                       - change),
                              "\n  ",
                              gettextf("Moved object %d from cluster %d to %d.",
                                       o, i, j)))
                fv_ids[o] <- j
                fv_p[i, ] <-
                    (fv_sums[i, ] - c(x[o, ])) / new_norms_rem[o]
                fv_p[j, ] <-
                    (fv_sums[j, ] + c(x[o, ])) / new_norms_add[o, j]
                fv_norms[c(i, j)] <-
                    c(new_norms_rem[o], new_norms_add[o, j])
                fv_sums[c(i, j), ] <-
                    fv_norms[c(i, j)] * fv_p[c(i, j), ]
                crossprods[, c(i, j)] <-
                    2 * g_tcrossprod(x, fv_sums[c(i, j), ])
                ## Record and update if optimal and big enough.
                if(change > opt_change) {
                    opt_change <- change
                    if(change > tol) {
                        opt_chains <- chains
                        ids <- fv_ids
                        sums <- fv_sums
                        norms <- fv_norms
                        p <- fv_p
                    }
                }
                chains <- chains + 1L
            }

            if(verbose && (opt_chains < maxchains))
                ## Indicate where the optimal chain occurred.
                message(gettextf("Rolling back %d chain move(s).",
                                 maxchains - opt_chains))

            ## If we found no change big enough, we are done.
            if(opt_chains == 0L)
                break

            ## Update value and resume iteration.
            new_value <- sum(norms)            
            old_value <- new_value
            iter <- iter + 1L
        }

        if(new_value > opt_value) {
            opt_value <- new_value
            opt_ids <- ids
            opt_p <- p
        }

        if(run >= nruns) break
        
        run <- run + 1L
        if(verbose)
            message(gettextf("Pclust run: %d", run))
        p <- start[[run]]
    }

    .hard_skmeans_object_for_normalized_x(x, opt_ids, k,
                                          p = opt_p, v = s - opt_value)
}

### * .skmeans_soft_pclust

.skmeans_soft_pclust <- 
function(x, k, m, weights = 1, control = NULL)
{
    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 100L
    
    reltol <- control$reltol
    if(is.null(reltol))
        reltol <- sqrt(.Machine$double.eps)
    
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")

    nr <- nrow(x)
    nc <- ncol(x)

    if(all(weights == 1)) {
        .W <- identity
    } else {
        if(!all(weights > 0))
            stop("All weights must be positive.")
        .W <- function(x) weights * x
    }

    ## Initialize.
    start <- control$start
    nruns <- length(start)
    
    opt_value <- Inf
    run <- 1L

    if(verbose && (nruns > 1L))
        message(gettextf("Pclust run: %d", run))

    repeat {
        p <- .skmeans_init_for_normalized_x(x, k, start[[run]], weights)
        old_value <- Inf
        iter <- 1L
        while(iter <= maxiter) {
            ## Fixed point iteration.
            dissimilarities <- pmax(1 - g_tcrossprod(x, p), 0)
            ## New memberhips.
            u <- clue:::.memberships_from_cross_dissimilarities(dissimilarities,
                                                                m)
            v <- .W(u ^ m)
            ## New prototypes.
            sums <- g_crossprod(v, x)
            norms <- row_norms(sums)
            p <- sums / norms
            ## New value.
            new_value <- sum(v) - sum(norms)
            
            if(verbose)
                message(gettextf("Iteration: %d *** value: %g",
                                 iter, new_value))

            if(abs(old_value - new_value)
                   < reltol * (abs(old_value) + reltol))
                    break
            old_value <- new_value
            iter <- iter + 1L
        }
        if(new_value < opt_value) {
            opt_value <- new_value
            opt_u <- u
            opt_p <- p
        }

        if(run >= nruns) break
        
        run <- run + 1L
        if(verbose)
            message(gettextf("Pclust run: %d", run))
        p <- start[[run]]
    }

    opt_ids <- max.col(opt_u)
    
    ## Ensure that opt_u is a stochastic matrix.
    opt_u <- pmax(opt_u, 0)
    opt_u <- opt_u / rowSums(opt_u)
    opt_u <- cl_membership(as.cl_membership(opt_u), k)

    dimnames(opt_p) <- list(seq_len(k), colnames(x))
    dimnames(opt_u) <- list(rownames(x), seq_len(k))
    names(opt_ids) <- rownames(x)

    out <- pclust_object(prototypes = opt_p,
                         membership = opt_u,
                         cluster = opt_ids,
                         family = skmeans_family,
                         m = m,
                         value = opt_value)
    class(out) <- unique(c("skmeans", class(out)))

    out
}

### * .skmeans_soft_pclust_via_clue_pclust

.skmeans_soft_pclust_via_clue_pclust <- 
function(x, k, m, weights = 1, control = NULL)
{
    ## Internally, we always keep x and prototypes normalized.
    .D_for_normalized_x_and_prototypes <-
        function(x, prototypes)
            pmax(1 - g_tcrossprod(x, prototypes), 0)
    .C_for_normalized_x_and_prototypes <-
        function(x, weights, control) {
            out <- g_col_sums_with_weights(x, weights)
            out / sqrt(sum(out ^ 2))
        }

    family <- skmeans_family
    C <- family$C
    D <- family$D
    family$C <- .C_for_normalized_x_and_prototypes
    family$D <- .D_for_normalized_x_and_prototypes
    out <- pclust(x, k, family, m, weights, control)
    out$family$C <- C
    out$family$D <- D

    class(out) <- unique(c("skmeans", class(out)))
    out
}
    
### * .skmeans_genetic

.skmeans_genetic <-
function(x, k, weights = NULL, control = NULL)
{
    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 12L
    
    mutations <- control$mutations
    if(is.null(mutations))
        mutations <- 0.1
    
    popsize <- control$popsize
    if(is.null(popsize))
        popsize <- 6L
    
    reltol <- control$reltol
    if(is.null(reltol))
        reltol <- sqrt(.Machine$double.eps)
    
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")

    nr <- nrow(x)
    nc <- ncol(x)

    ## Normalize x.
    x <- row_normalize(x)

    ## In the weighted case, we can perform all computations on w_i x_i
    ## (after normalization).
    if(all(weights == 1)) {
        s <- nr
    } else {
        if(!all(weights > 0))
            stop("All weights must be positive.")
        s <- sum(weights)
        x <- weights * x
    }

    ## Initialize p.
    start <- control$start
    if(is.null(start))
        start <- as.list(rep.int("p", popsize))
    else if(!is.list(start) || inherits(start, "skmeans"))
        start <- list(start)
    p <- lapply(start,
                function(s)
                .skmeans_init_for_normalized_x(x, k, s))
    popsize <- length(p)
    
    ## Initialize ids.
    ids <- lapply(p,
                  function(p1)
                  ids_from_similarities(g_tcrossprod(x, p1), k))

    if(verbose)
        message(gettextf("Initial solutions created (length: %d)",
                         popsize))

    value <- rep.int(0, popsize)

    genetic_mutator <- function(id) {
        ## Randomly change cluster membership.
        ## Simplistic solution.
        mutations2 <- mutations * k / (k - 1L)
        newid <- id
        mut <- ceiling(runif(nr * mutations2, 0, nr))
        newid[mut] <- ceiling(runif(nr * mutations2, 0, k))
        newid
    }

    genetic_selection <- function(v) {
        v <- v + runif(length(v), min(v), max(v))
        which(v >= sort(v, decreasing = TRUE)[popsize])
    }

    iter <- 1L
    while(iter <= maxiter) {
        ## Select best prototypes.
        sel <- genetic_selection(value)
        if(verbose)
            message(gettextf("Selecting new population: %s",
                             paste(sel, collapse = ",")))
        ids <- ids[sel]
        p <- p[sel]
        value <- value[sel]
        for (i in seq_len(popsize)) {
            ## Generate new prototypes by mutating and k-means.
            new <- popsize + i
            ## <NOTE>
            ids[[new]] <- genetic_mutator(ids[[i]])
            ## Currently, genetic_mutator() does not ensure that all
            ## class ids are used: .hard_skmeans_C_for_normalized_x()
            ## will add unused ones back if needed.
            sums <- .hard_skmeans_C_for_normalized_x(x, ids[[new]])
            ## </NOTE>
            norms <- row_norms(sums)
            if(verbose)
                message(gettextf("Iteration: %d *** value[%d]: %g pre-optimized",
                                 iter, new, s - sum(norms)))
            repeat {
                p[[new]] <- sums / norms
                similarities <- g_tcrossprod(x, p[[new]])
                ids[[new]] <- ids_from_similarities(similarities, k)
                sums <- .hard_skmeans_C_for_normalized_x(x, ids[[new]])
                norms <- row_norms(sums)
                oldvalue <- value[new]
                value[new] <- sum(norms)
                if(!is.na(oldvalue) &&
                    abs(oldvalue - value[new]) < reltol*(oldvalue + reltol)) break
            }
            if(verbose)
                message(gettextf("Iteration: %d *** value[%d]: %g",
                                 iter, new, s - value[new]))

        }
        iter <- iter + 1L
    }

    ## Get ids from the winner population.
    ids <- ids[[which.max(value)]]

    .hard_skmeans_object_for_normalized_x(x, ids, k, s)
}

### * .skmeans_local_improvement_heuristic

## <FIXME>
## Remove eventually.
.skmeans_local_improvement_heuristic <-
function(x, k, control = NULL)
{
    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 100L
    reltol <- control$reltol
    if(is.null(reltol))
        reltol <- 1e-8
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")

    nr <- nrow(x)
    nc <- ncol(x)

    ## Normalize x.
    x <- row_normalize(x)

    ## Initialize p.
    p <- .hard_skmeans_init_for_normalized_x(x, k, control)

    old_value <- 0

    iter <- 1L
    while(iter <= maxiter) {
        similarities <- g_tcrossprod(x, p)
        ids <- ids_from_similarities(similarities, k)
        ## New prototypes.
        sums <- .hard_skmeans_C_for_normalized_x(x, ids)
        norms <- row_norms(sums)
        p <- sums / norms
        new_value <- sum(norms)
        if(verbose)
            message(gettextf("Iteration: %d [KM] *** value: %g",
                             iter, nr - new_value))
        if(abs(old_value - new_value)
           < reltol * (abs(old_value) + reltol)) {
            ## Block update performed only minor improvement.
            ## Try first variation steps.
            crossprods <- 2 * sweep(similarities, 2L, norms, "*")
            count <- 1L
            repeat {
                nids <- norms[ids]
                ind <- cbind(seq_len(nr), ids)
                ## The effects of removing an object from its cluster.
                new_norms_rem <- sqrt(nids ^ 2 - crossprods[ind] + 1)
                Delta_Q_rem <- new_norms_rem - nids
                ## The effects of adding an object to another cluster.
                new_norms_add <-
                    sqrt(sweep(crossprods, 2, norms ^ 2 + 1, "+"))
                Delta_Q_add <- sweep(new_norms_add, 2, norms, "-")
                ## What is the best such move?
                Delta_Q <- Delta_Q_rem + Delta_Q_add
                Delta_Q[ind] <- 0
                Delta_Q <- as.numeric(Delta_Q)
                pos <- which.max(Delta_Q)
                ## This changes object o from its cluster i (ids[o]) to
                ## cluster j.  If the change is large enough, perform
                ## it; otherwise we are done.
                if(Delta_Q[pos] < reltol * (abs(new_value) + reltol))
                    break
                j <- (pos - 1L) %/% nr + 1L
                o <- (pos - 1L) %% nr + 1L
                i <- ids[o]
                ids[o] <- j
                if(verbose)
                    message(gettextf("Moving %d from cluster %d to %d.",
                                     o, i, j))
                p[i, ] <- (sums[i, ] - c(x[o, ])) / new_norms_rem[o]
                p[j, ] <- (sums[j, ] + c(x[o, ])) / new_norms_add[o, j]
                new_value <- new_value + Delta_Q[pos]
                if(verbose)
                    message(gettextf("Iteration: %d [FV %d] *** value: %g",
                                     iter, count, nr - new_value))
                count <- count + 1L
                ## <NOTE>
                ## Of course, this is rather inefficient, but let's just
                ## see if we can find better solutions by repeating
                ## first variation ...
                norms[c(i, j)] <-
                    c(new_norms_rem[o], new_norms_add[o, j])
                sums[c(i, j), ] <-
                    norms[c(i, j)] * p[c(i, j), ]
                crossprods[, c(i, j)] <-
                    2 * g_tcrossprod(x, sums[c(i, j), ])
                ## </NOTE>
            }
            if(count == 1L) break
        }
        old_value <- new_value
        iter <- iter + 1L
    }

    .hard_skmeans_object_for_normalized_x(x, ids, k, nr)
}
## </FIXME>

### * .skmeans_local_improvement_heuristic_with_chains

## <FIXME>
## Remove eventually.
.skmeans_local_improvement_heuristic_with_chains <-
function(x, k, control = NULL)
{
    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 100L
    maxchains <- control$maxchains
    if(is.null(maxchains))
        maxchains <- 10L
    reltol <- control$reltol
    if(is.null(reltol))
        reltol <- 1e-8
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")

    nr <- nrow(x)
    nc <- ncol(x)

    ## Normalize x.
    x <- row_normalize(x)

    ## Initialize p.
    p <- .hard_skmeans_init_for_normalized_x(x, k, control)

    old_value <- 0

    iter <- 1L
    while(iter <= maxiter) {
        similarities <- g_tcrossprod(x, p)
        ids <- ids_from_similarities(similarities, k)
        ## New prototypes.
        sums <- .hard_skmeans_C_for_normalized_x(x, ids)
        norms <- row_norms(sums)
        p <- sums / norms
        new_value <- sum(norms)
        if(verbose)
            message(gettextf("Iteration: %d [KM] *** value: %g",
                             iter, nr - new_value))
        if(abs(old_value - new_value)
           < reltol * (abs(old_value) + reltol)) {
            ## Block update performed only minor improvement.
            ## Try first variation steps.
            crossprods <- 2 * sweep(similarities, 2L, norms, "*")
            count <- 1L
            repeat {
                ## Save prototypes for rollback.
                oldp <- p
                change <- 0L
                marked <- integer()
                chains <- 1L
                ## Try maxchains steps before discarding first variation
                ## (Kernighan-Lin chains)
                while (chains <= maxchains) {
                    nids <- norms[ids]
                    ind <- cbind(seq_len(nr), ids)
                    ## The effects of removing an object from its
                    ## cluster.
                    new_norms_rem <-
                        sqrt(nids ^ 2 - crossprods[ind] + 1)
                    Delta_Q_rem <- new_norms_rem - nids
                    ## The effects of adding an object to another
                    ## cluster.
                    new_norms_add <-
                        sqrt(sweep(crossprods, 2L, norms ^ 2 + 1, "+"))
                    Delta_Q_add <- sweep(new_norms_add, 2, norms, "-")
                    ## What is the best such move?
                    Delta_Q <- Delta_Q_rem + Delta_Q_add
                    Delta_Q[ind] <- -1e9
                    Delta_Q[marked,] <- -1e9
                    Delta_Q <- as.numeric(Delta_Q)
                    pos <- which.max(Delta_Q)
                    ## This changes object o from its cluster i (ids[o])
                    ## to cluster j.  If the change is large enough,
                    ## perform it; otherwise we are done.
                    j <- (pos - 1L) %/% nr + 1L
                    o <- (pos - 1L) %% nr + 1L
                    marked <- c(marked, o)
                    i <- ids[o]
                    ids[o] <- j
                    if(verbose)
                        message(gettextf("Moving %d from cluster %d to %d.",
                                         o, i, j))
                    p[i, ] <- (sums[i, ] - c(x[o, ])) / new_norms_rem[o]
                    p[j, ] <- (sums[j, ] + c(x[o, ])) / new_norms_add[o, j]
                    change <- Delta_Q[pos] + change
                    if(verbose)
                        message(gettextf("Iteration: %d [FV %d, CH %d] *** value: %g change: %g",
                                         iter, count, chains,
                                         nr - new_value, - 2 * change))

                    ## <NOTE>
                    ## Of course, this is rather inefficient, but let's
                    ## just see if we can find better solutions by
                    ## repeating first variation ...
                    norms[c(i, j)] <-
                        c(new_norms_rem[o], new_norms_add[o, j])
                    sums[c(i, j), ] <-
                        norms[c(i, j)] * p[c(i, j), ]
                    crossprods[, c(i, j)] <-
                        2 * g_tcrossprod(x, sums[c(i, j), ])
                    ## </NOTE>
                    if(change > reltol * (abs(new_value) + reltol)) {
                        count <- count + 1L
                        break
                    }
                    chains <- chains + 1L
                }
                if (chains == maxchains + 1L) {
                    p <- oldp
                    if (verbose)
                        message(gettextf("Rolling back %d chain moves.",
                                         chains - 1L))
                    break
                }
                new_value <- sum(norms)
            }
            if(count == 1L) break
        }
        old_value <- new_value
        iter <- iter + 1L
    }

    ## Fix ids from chain evaluation.
    ids <- ids_from_similarities(g_tcrossprod(x, p), k)

    .hard_skmeans_object_for_normalized_x(x, ids, k, nr)
}
## </FIXME>

### * .skmeans_CLUTO

## Spherical sparse k-means via CLUTO file I/O.

.skmeans_CLUTO <-
function(x, k, control = NULL)
{
    ## Try to determine path to CLUTO vcluster executable.
    vcluster <- control$vcluster
    if(is.null(vcluster))
        vcluster <- "vcluster"
    if(!file.exists(Sys.which(vcluster)))
        stop("CLUTO vcluster executable not found")

    colmodel <- control$colmodel
    if(is.null(colmodel))
        colmodel <- "none"
    
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")
    
    ifile <- control$ifile

    ## Could add more fancy control list expansion eventually.
    control <- paste(as.character(control$control), collapse = " ")

    tmp <- tempfile()
    if (is.null(ifile)) {
        datfile <- tmp
        on.exit(file.remove(datfile))
        write_stm_CLUTO(x, datfile)
    }
    else
        datfile <- ifile

    clustfile <- paste(tmp, ".clustering.", k, sep = "")
    on.exit(file.remove(clustfile), add = TRUE)

    cmdout <-
        system(sprintf("%s -colmodel=%s -clustfile=%s %s %s %d",
                       vcluster, colmodel, clustfile, control,
                       datfile, k),
               intern = TRUE)
    if(verbose)
        message(paste(cmdout, collapse = "\n"))

    result <- read.table(clustfile, header = FALSE) + 1
    ids <- result$V1

    .hard_skmeans_object_for_normalized_x(row_normalize(x), ids, k)
}

### * .skmeans_gmeans

.skmeans_gmeans <-
function(x, k, control = NULL)
{
    gmeans <- control$gmeans
    if(is.null(gmeans))
        gmeans <- "gmeans"
    if(!file.exists(Sys.which(gmeans)))
        stop("gmeans executable not found")

    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")

    ifile <- control$ifile

    ## Vector of ids
    start <- control$start
    if (!is.null(start)) {
        initfile <- tempfile()
        on.exit(file.remove(initfile))
        write(c(length(start), start - 1), initfile, sep = "\n")
    }

    ## Could add more fancy control list expansion eventually.
    control <- paste(as.character(control$control), collapse = " ")

    tmp <- tempfile()
    if (is.null(ifile)) {
        datfile <- tmp
        on.exit(file.remove(sprintf("%s_dim", datfile),
                            sprintf("%s_col_ccs", datfile),
                            sprintf("%s_row_ccs", datfile),
                            sprintf("%s_tfn_nz", datfile)),
                add = TRUE)
        slam::write_stm_MC(x, datfile)
    }
    else
        datfile <- ifile

    clustfile <- paste(tmp, "_tfn_doctoclus.", k, sep = "")
    on.exit(file.remove(clustfile), add = TRUE)

    cmdout <-
        system(sprintf("%s -c %s %s -O %s %s %s",
                       gmeans, k,
                       if(!is.null(start)) paste("-i f", initfile) else "",
                       tmp, control, datfile),
               intern = TRUE)
    if(verbose)
        message(paste(cmdout, collapse = "\n"))

    result <- read.table(clustfile) + 1
    ## The number in the first row is the number of data items clustered
    ids <- result$V1[-1L]

    .hard_skmeans_object_for_normalized_x(row_normalize(x), ids, k)
}

### .skmeans_kmndirs

.skmeans_kmndirs <-
function(x, k, control = NULL)
{
    nrandom <- control$nrandom
    if(is.null(nrandom))
        nrandom <- 1000L

    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 10L

    x <- row_normalize(x)

    ids <- kmndirs:::kmndirs(x, k, nrandom, maxiter) + 1

    .hard_skmeans_object_for_normalized_x(x, ids, k)
}

## * Helpers

## In the "hard" case (m = 1) with normalized x:
## * the consensus function is the weighted sum of all objects in a
##   group;
## * the value is the sum of the weights minus the sum of the norms of
##   the consensus sums.
## Note that we always replace x by weights * x in the R-based fitters.

## Vectorized consensus for the hard normalized x case:

.hard_skmeans_C_for_normalized_x <-
function(x, ids)
{
    ## Somewhat costly ...
    ids <- factor(ids)
    
    out <- g_col_sums_by_group(x, ids)

    ## An id of 0 actually indicates missingness (from CLUTO),
    ## so needs to be dropped.
    if(!is.na(pos <- match(0, levels(ids))))
        out <- out[-pos, , drop = FALSE]

    out
}

## Initializer for the hard normalized x case.

## <FIXME>
## Remove eventually.
.hard_skmeans_init_for_normalized_x <-
function(x, k, control)
{
    p <- control$start
    if(inherits(p, "skmeans")) {
        p <- p$prototypes
    } else {
        ## In case this is a list:
        if(is.list(p))
            p <- p[[1L]]
        ## In case this got us nowhere:
        if(is.null(p)) {
            ## Initialize p by choosing k random rows of x.
            ## (Hence initial prototypes are already normalized.)
            p <- as.matrix(x[sample.int(nrow(x), k), , drop = FALSE])
        }
    }
    p
}
## </FIXME>

## Initializer for normalized x.

.skmeans_init_for_normalized_x <-
function(x, k, start, weights = 1)
{
    if(is.character(start)) {
        if(start == "p") {
            ## Initialize prototypes by choosing k random rows of x.
            ## (Hence prototypes are already normalized.)
            as.matrix(x[sample.int(nrow(x), k), , drop = FALSE])
        }
        else if(start == "i") {
            ## Initialize by choosing random ids.
            ids <- sample.int(k, nrow(x), replace = TRUE)
            ## Could ensure that all ids are used.
            .hard_skmeans_C_for_normalized_x(x, ids)
        }
        else if(start == "S") {
            p1 <- row_normalize(rbind(g_col_sums(weights * x)))
            .skmeans_init_helper(x, k, p1)
        }
        else if(start == "s") {
            p1 <- x[sample.int(nrow(x), 1L), , drop = FALSE]
            .skmeans_init_helper(x, k, p1)
        }
        else
            stop(gettextf("Invalid control option 'start'"),
                 domain = NA)
    }
    else if(inherits(start, "skmeans"))
        row_normalize(start$prototypes)
    else if(!is.null(dim(start)))
        row_normalize(start)
    else {
        ## A vector of class ids, hopefully.
        ids <- match(start, unique(start))
        .hard_skmeans_C_for_normalized_x(x, ids)
    }
}

.skmeans_init_helper <-
function(x, k, p1)
{
    ## For normalized x and first prototype p1, determine k - 1 further
    ## prototypes as "items" which are farthest away (least similar) to
    ## all previously picked prototypes.
    p <- matrix(0, k, ncol(x))
    p1 <- rbind(c(as.matrix(p1)))
    s <- g_tcrossprod(x, p1)
    p[1L, ] <- p1
    for(i in seq(2L, length.out = k - 1L)) {
        s <- pmax(s, g_tcrossprod(x, p1))
        p[i, ] <- p1 <- as.matrix(x[which.min(s), , drop = FALSE])
    }
    p
}

## Object generator for the hard normalized x case.

.hard_skmeans_object_for_normalized_x <-
function(x, ids, k, s = nrow(x), p = NULL, v = NULL)
{
    if(is.null(p) || is.null(v)) {
        ## If we only have the normalized prototypes we cannot figure
        ## out the value of the criterion function.
        sums <- .hard_skmeans_C_for_normalized_x(x, ids)
        norms <- row_norms(sums)
        p <- sums / norms
        v <- s - sum(norms)
    }

    ## Handle 0 ids indicating missingness (CLUTO).
    if(any(ind <- (ids == 0))) {
        ids[ind] <- max.col(g_tcrossprod(x[ind, , drop = FALSE], p))
    }

    ## <NOTE>
    ## No longer provide redundant memberships for hard partitions.
    ## </NOTE>

    k <- nrow(p)                        # Just making sure ...
    dimnames(p) <- list(seq_len(k), colnames(x))
    names(ids) <- rownames(x)

    out <- pclust_object(prototypes = p,
                         membership = NULL,
                         cluster = ids,
                         family = skmeans_family,
                         m = 1,
                         value = v)
    class(out) <- unique(c("skmeans", class(out)))
    
    out
}

## I2 for the hard case.

I2 <-
function(x, ids)
{
    ids <- match(ids, unique(ids))
    sums <- .hard_skmeans_C_for_normalized_x(row_normalize(x), ids)
    sum(row_norms(sums))
}

### * Utilities

row_norms <-
function(x)
    sqrt(g_row_sums(x ^ 2))

row_normalize <-
function(x)
    x / row_norms(x)

g_col_sums_with_weights <-
function(x, w)
{
    ## Improve performance by leaving out Ops dispatch.
    ## (Other sparse matrix classes could be dealt with similarly ...)
    if(inherits(x, "simple_triplet_matrix")) {
        x$v <- x$v * w[x$i]
        slam::col_sums(x)
    } else if(inherits(x, "dgCMatrix")) {
        x@x <- x@x * w[x@i + 1L]
        Matrix::colSums(x)
    } else {
        g_col_sums(w * x)
    }
}

g_col_sums_with_logical_index <-
function(x, l)
{
    ## Improve performance by leaving out Ops dispatch.
    ## (Other sparse matrix classes could be dealt with similarly ...)
    if(inherits(x, "simple_triplet_matrix")) {
        x$v <- x$v * l[x$i]
        slam::col_sums(x)
    } else if(inherits(x, "dgCMatrix")) {
        x@x <- x@x * l[x@i + 1L]
        Matrix::colSums(x)
    } else {
        g_col_sums(x[l, , drop = FALSE])
    }
}

ids_from_similarities <-
function(x, k)
{
    ids <- max.col(x)
    all_ids_used <- sort(unique(ids))
    if(length(all_ids_used) < k) {
        ## Assign objects with the smallest similarities to their own
        ## cluster.
        o <- order(x[cbind(seq_along(ids), ids)])
        unused <- setdiff(seq_len(k), all_ids_used)
        ids[o[seq_along(unused)]] <- unused
    }
    ids
}


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
