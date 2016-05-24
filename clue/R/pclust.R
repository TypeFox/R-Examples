### * cl_pclust

cl_pclust <-
function(x, k, method = NULL, m = 1, weights = 1, control = list())
{
    ## Partition a cluster ensemble x into (at most) k classes by
    ## minimizing
    ##   \sum_b \sum_j w_b u_{bj}^m d(x_b, p_j) ^ e
    ## for "suitable" prototypes p_1, ..., p_k, where 1 <= m < \infty,
    ## with 1 corresponding to hard (secondary) partitions, and d a
    ## dissimilarity measure (such as Euclidean dissimilarity of
    ## partitions or hierarchies).
    ##
    ## The algorithm works whenever there is a consensus method for
    ## solving
    ##   \sum_b u_{bj}^m d(x_b, p) ^ e => \min_p
    ##
    ## As we refer to consensus methods by their *name* (e.g., 'HBH'),
    ## we rely on the registration mechanism (set_cl_consensus_method())
    ## to provide the required information about d and e.

    clusterings <- as.cl_ensemble(x)

    type <- .cl_ensemble_type(clusterings)

    if(type == "partition") {
        ## Canonicalize by turning into an ensemble of partitions
        ## represented by membership matrices with the same (minimal)
        ## number of columns.
        memberships <-
            lapply(clusterings, cl_membership,
                   max(sapply(clusterings, n_of_classes)))
        clusterings <-
            cl_ensemble(list = lapply(memberships, as.cl_partition))
    }

    if(!inherits(method, "cl_consensus_method")) {
        ## Get required information on d and e from the registry.
        if(is.null(method))
            method <- .cl_consensus_method_default(type)
        method <- get_cl_consensus_method(method, type)
        ## Note that this avoids registry lookup in subsequent calls to
        ## cl_consensus().
        if(is.null(method$exponent))
            stop("No information on exponent in consensus method used.")
        e <- method$exponent
        if(is.null(method$dissimilarity))
            stop("No information on dissimilarity in consensus method used.")
        d <- function(x, y = NULL)
            cl_dissimilarity(x, y, method = method$dissimilarity)
    }

    family <- pclust_family(D = d, C = method$definition, e = e)
    out <- pclust(x, k, family, m, weights, control)

    ## Massage the results a bit.
    dissimilarities <- as.matrix(d(clusterings) ^ e)
    out$call <- match.call()
    out <- .structure(c(out,
                        list(silhouette =
                             silhouette(out$cluster,
                                        dmatrix = dissimilarities),
                             validity =
                             cl_validity(cl_membership(out),
                                         dissimilarities),
                             ## <NOTE>
                             ## Information about d and e is also in the
                             ## family returned, of course.  Trying to be
                             ## nice to users by "directly" providing d
                             ## and e is currently of limited usefulness
                             ## as the pclust representation is not
                             ## directly available to users.
                             d = d,
                             e = e
                             ## </NOTE>
                             )),
                      class = unique(c("cl_pclust", class(out))))

    as.cl_partition(out)
}

print.cl_pclust <-
function(x, ...)
{
    txt <- if(x$m == 1)
        gettextf("A hard partition of a cluster ensemble with %d elements into %d classes.",
                 n_of_objects(x), n_of_classes(x))
    else
        gettextf("A soft partition (degree m = %f) of a cluster ensemble with %d elements into %d classes.",
                x$m, n_of_objects(x), n_of_classes(x))
    writeLines(strwrap(txt))
    NextMethod("print", x, header = FALSE)
    print(x$validity, ...)
    invisible(x)
}

### * pclust

pclust <-
function(x, k, family, m = 1, weights = 1, control = list())
{
    ## A general purpose alternating optimization algorithm for
    ## prototype-based partitioning.
    
    ## For now, assume family specifies three functions:
    ## * A dissimilarity function D() for data and prototypes.
    ## * A consensus function C() for data, weights and control.
    ## * An init function init() of data and k giving an initial object
    ##   of k prototypes.
    ##
    ## We use k as the second argument as this seems to be common
    ## practice for partitioning algorithms.

    ## <NOTE>
    ## We assume that consensus functions can all handle WEIGHTS
    ## (formals: x, weights, control; only used positionally).
    ## <NOTE>
    
    ## <NOTE>
    ## We now allow for arbitrary representations/objects of prototypes.
    ## What is needed are functions to modify a *single* prototype and
    ## subset the prototypes.  By default, list and matrix (with the
    ## usual convention that rows are "objects") representations are
    ## supported.  Otherwise, the family needs to provide suitable
    ## .modify() and .subset() functions.
    ## The approach relies on having the initializer of the family
    ## (init()) return an appropriate object of prototypes.
    ## It would be possible to have default initializers as well to
    ## randomly subset the data (i.e., select elements of lists or rows
    ## of matrices, respectively).
    ## </NOTE>

    ## <NOTE>
    ## The 'prototypes' are not necessarily objects of the same kind as
    ## the data objects.  Therefore, D() is really a 2-argument
    ## cross-dissimilarity function.
    ## It would also be useful to have a way of computing the pairwise
    ## dissimilarities between objects: but this is something different
    ## from D() is objects and prototypes are not of the same kind.
    ## A "clean" solution could consist in specifying the family either
    ## via a (non-symmetric) cross-dissimilarity function X(), or a
    ## symmetric D() which when called with a single argument gives the
    ## pairwise object dissimilarities.
    ## I.e., 
    ##   pclust_family(D = NULL, C, init = NULL, X = NULL, ......)
    ## using
    ## * If both D and X are not given => TROUBLE.
    ## * If only D is given: use for X as well.
    ## * If only X is given: only use as such.
    ## Something for the future ...
    ## </NOTE>

    ## <NOTE>
    ## If people have code for computing cross-dissimilarities for the
    ## data and a *single* prototype (say, xd()), they can easily wrap
    ## into what is needed using
    ##   t(sapply(prototypes, function(p) xd(x, p)))
    ## Assuming symmetry of the dissimilarity, they could also do
    ##   t(sapply(prototypes, xd, x))
    ## </NOTE>

    ## Perhaps check whether 'family' is a feasible/suitable pclust
    ## family (object).
    D <- family$D
    C <- family$C
    e <- family$e
    .modify <- family$.modify
    .subset <- family$.subset

    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 100L
    nruns <- control$nruns
    reltol <- control$reltol
    if(is.null(reltol))
        reltol <- sqrt(.Machine$double.eps)
    start <- control$start
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")
    ## Do this at last ...
    control <- as.list(control$control)

    ## Handle start values and number of runs.
    if(!is.null(start)) {
        if(!is.list(start)) {
            ## Be nice to users.
            start <- list(start)
        }
        nruns <- length(start)
    } else {
        if(is.null(nruns)) {
            ## Use nruns only if start is not given.
            nruns <- 1L
        }
        start <- replicate(nruns, family$init(x, k), simplify = FALSE)
    }
    
    ## Initialize.
    ## We need to do this here because it is (currently) the only way to
    ## figure out the number B of objects to be partitioned (which is
    ## needed for getting the object weights to the right length).
    prototypes <- start[[1L]]
    dissimilarities <- D(x, prototypes) ^ e
    B <- NROW(dissimilarities)
    ## Also try to figure out (if necessary) how to modify a single
    ## prototype and to subset the prototypes.  Note that we can only
    ## check this *after* prototypes were obtained (and not when the
    ## family object is created).
    if(is.null(.modify)) {
        if(is.list(prototypes))
            .modify <- function(x, i, value) {
                x[[i]] <- value
                x
            }
        else if(is.matrix(prototypes))
            .modify <- function(x, i, value) {
                x[i, ] <- value
                x
            }
        else stop("Cannot determine how to modify prototypes.")
    } else if(!is.function(.modify) ||
            !identical(formals(args(.modify)), c("x", "i", "value")))
        stop("Invalid function to modify prototypes.")
    if(is.null(.subset)) {
        if(is.list(prototypes))
            .subset <- `[`
        else if(is.matrix(prototypes))                 
            .subset <- function(x, i)
                x[i, , drop = FALSE]
        else stop("Cannot determine how to subset prototypes.")
    } else if(!is.function(.subset) ||
            !identical(formals(args(.subset)), c("x", "i")))
        stop("Invalid function to subset prototypes.")

    weights <- rep(weights, length.out = B)
    if(any(weights < 0))
        stop("Argument 'weights' has negative elements.")
    if(!any(weights > 0))
        stop("Argument 'weights' has no positive elements.")

    ## A little helper.
    .make_unit_weights <- function(B, i) {
        out <- double(B)
        out[i] <- 1
        out
    }

    if(m == 1) {
        ## Hard partitions.
        value <- if(all(weights == 1))
            function(dissimilarities, ids)
                sum(.one_entry_per_column(dissimilarities, ids))
        else
            function(dissimilarities, ids)
                sum(weights *
                    .one_entry_per_column(dissimilarities, ids))
        opt_value <- Inf
        run <- 1L
        if(verbose && (nruns > 1L))
            message(gettextf("Pclust run: %d", run))
        repeat {
            class_ids <- max.col( - dissimilarities )
            old_value <- value(dissimilarities, class_ids)
            if(verbose)
                message(gettextf("Iteration: 0 *** value: %g", old_value))
            iter <- 1L
            while(iter <= maxiter) {
                class_ids_used <- unique(class_ids)
                for(j in class_ids_used)
                    prototypes <-
                        .modify(prototypes, j, 
                                C(x, weights * (class_ids %in% j),
                                  control))
                dissimilarities <- D(x, prototypes) ^ e
                class_ids <- max.col( - dissimilarities )
                ## Try avoiding degenerate solutions.
                if(length(class_ids_used) < k) {
                    ## Find the k - l largest
                    ## object-to-assigned-prototype dissimilarities.
                    o <- order(.one_entry_per_column(dissimilarities,
                                                     class_ids),
                               decreasing = TRUE)
                    ## Find and recompute unused prototypes.
                    unused <- setdiff(seq_len(k), class_ids_used)
                    for(j in seq_along(unused))
                        prototypes <-
                            .modify(prototypes, unused[j],
                                    C(x, .make_unit_weights(B, o[j]),
                                      control))
                    dissimilarities[, unused] <-
                        D(x, .subset(prototypes, unused)) ^ e
                    class_ids <- max.col( - dissimilarities )
                    ## For the time being, do not retry in case the
                    ## solution is still degenerate.
                }
                new_value <- value(dissimilarities, class_ids)
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
                converged <- (iter <= maxiter)
                opt_value <- new_value
                opt_class_ids <- class_ids
                opt_prototypes <- prototypes
            }
            if(run >= nruns) break
            run <- run + 1L
            if(verbose)
                message(gettextf("Pclust run: %d", run))
            prototypes <- start[[run]]
            dissimilarities <- D(x, prototypes) ^ e
        }
        ## We should really have a suitable "sparse matrix" class for
        ## representing the memberships of hard partitions.  For now:
        opt_u <- NULL
        ## opt_u <- matrix(0, B, k)
        ## opt_u[cbind(seq_len(B), opt_class_ids)] <- 1
    }
    else {
        ## Soft partitions.
        value <- if(all(weights == 1))
            function(dissimilarities, u)
                sum(u ^ m * dissimilarities)
        else
            function(dissimilarities, u)
                sum(weights *
                    u ^ m * dissimilarities)
        opt_value <- Inf
        run <- 1L
        if(verbose && (nruns > 1L))
            message(gettextf("Pclust run: %d", run))
        repeat {
            u <- .memberships_from_cross_dissimilarities(dissimilarities,
                                                         m)
            old_value <- value(dissimilarities, u)
            if(verbose)
                message(gettextf("Iteration: 0 *** value: %g", old_value))
            iter <- 1L
            while(iter <= maxiter) {
                ## Update the prototypes.
                ## This amounts to solving, for each j:
                ##   \sum_b w_b u_{bj}^m D(x_b, p) ^ e => \min_p
                ## I.e., p_j is the *weighted* consensus of the x_b with
                ## corresponding weights u_{bj}^m.
                for(j in seq_len(k)) {
                    prototypes <-
                        .modify(prototypes, j,
                                C(x, weights * u[, j] ^ m, control))
                }
                ## Update u.
                dissimilarities <- D(x, prototypes) ^ e
                u <- .memberships_from_cross_dissimilarities(dissimilarities,
                                                             m)
                new_value <- value(dissimilarities, u)
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
                converged <- (iter <= maxiter)                
                opt_value <- new_value
                opt_prototypes <- prototypes
                opt_u <- u
            }
            if(run >= nruns) break
            run <- run + 1L
            if(verbose)
                message(gettextf("Pclust run: %d", run))
            prototypes <- start[[run]]
            dissimilarities <- D(x, prototypes) ^ e
        }
        opt_class_ids <- max.col(opt_u)
        ## Ensure that opt_u is a stochastic matrix.
        opt_u <- pmax(opt_u, 0)
        opt_u <- opt_u / rowSums(opt_u)
        rownames(opt_u) <- rownames(dissimilarities)
        opt_u <- cl_membership(as.cl_membership(opt_u), k)
    }

    names(opt_class_ids) <- rownames(dissimilarities)

    pclust_object(prototypes = opt_prototypes,
                  membership = opt_u,
                  cluster = opt_class_ids,
                  family = family,
                  m = m,
                  value = opt_value,
                  call = match.call(),
                  attributes = list("converged" = converged))
}

print.pclust <-
function(x, header = TRUE, ...)
{
    is_hard <- (x$m == 1)
    class_ids <- cl_class_ids(x)
    if(header) {
        txt <- if(is_hard)
            gettextf("A hard partition of %d objects into %d classes.",
                     length(class_ids), n_of_classes(x))
        else
            gettextf("A soft partition (degree m = %f) of %d objects into %d classes.",
                     x$m, length(class_ids), n_of_classes(x))
        writeLines(strwrap(txt))
    }
    if(is_hard) {
        print(class_ids, ...)
    }
    else {
        writeLines("Class memberships:")
        print(cl_membership(x), ...)
        writeLines("Class ids of closest hard partition:")
        print(unclass(class_ids), ...)
    }
    invisible(x)
}
    

### * pclust_family

pclust_family <-
function(D, C, init = NULL, description = NULL, e = 1,
         .modify = NULL, .subset = NULL)
{
    ## Add checking formals (lengths) eventually ...
    if(is.null(init)) {
        ## Works for list representations ...
        init <- function(x, k) sample(x, k)
    }
    .structure(list(description = description,
                    D = D, C = C, init = init, e = e,
                    .modify = .modify, .subset = .subset),
               class = "pclust_family")
}

### * pclust_object

pclust_object <-
function(prototypes, membership, cluster, family, m = 1, value, ...,
         classes = NULL, attributes = NULL)
{
    out <- c(list(prototypes = prototypes,
                  membership = membership,
                  cluster = cluster,
                  family = family,
                  m = m,
                  value = value),
             list(...))
    attributes(out) <- c(attributes(out), attributes)
    classes <- unique(as.character(classes))
    class(out) <- c(classes[classes != "pclust"], "pclust")
    out
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
