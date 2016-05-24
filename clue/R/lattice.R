cl_meet <-
function(x, y)
{
    ## General case.
    ## x either an ensemble, or x and y two clusterings with the same
    ## number of objects.
    if(!inherits(x, "cl_ensemble")) {
        ## Be nice about error messages.
        if(n_of_objects(x) != n_of_objects(y))
            stop("Arguments 'x' and 'y' must have the same number of objects.")
        x <- cl_ensemble(x, y)
    }

    if(inherits(x, "cl_partition_ensemble"))
        .cl_meet_partition(x)
    else if(inherits(x, "cl_dendrogram_ensemble"))
        .cl_meet_dendrogram(x)
    else if(inherits(x, "cl_hierarchy_ensemble"))
        .cl_meet_hierarchy(x)
    else
        stop("Cannot compute meet of given clusterings.")
}

.cl_meet_partition <-
function(x)    
{
    x <- unique(x)
    if(length(x) == 1L)
        return(cl_partition_by_class_ids(cl_class_ids(x[[1L]])))

    ids <- seq_len(n_of_objects(x[[1L]]))
    ## Cross-classify the objects.
    z <- split(ids, lapply(x, cl_class_ids))
    ## Subscript on the non-empty cells to get adjacent class ids.
    lens <- sapply(z, length)
    pos <- which(lens > 0)
    ids[unlist(z, use.names = FALSE)] <-
        rep(seq_along(z[pos]), lens[pos])
    cl_partition_by_class_ids(ids)
}

.cl_meet_dendrogram <-
function(x)
{
    ## Meet of an ensemble of dendrograms.
    ## We need the maximal ultrametric dominated by the given ones,
    ## which can be obtained by hierarchical clustering with single
    ## linkage on the pointwise minima of the ultrametrics.
    as.cl_dendrogram(hclust(as.dist(do.call(pmin,
                                            lapply(x, cl_ultrametric))),
                            "single"))
}

.cl_meet_hierarchy <-
function(x)
{
    ## Meet of an ensemble of n-trees.
    ## Need to find the classes in *all* n-trees.
    ## Equivalent to computing a strict majority tree.
    .cl_consensus_hierarchy_majority(x, rep.int(1, length(x)),
                                     list(p = 1))
}

cl_join <-
function(x, y)
{
    ## General case.
    ## x either an ensemble, or x and y two clusterings with the same
    ## number of objects.
    if(!inherits(x, "cl_ensemble")) {
        ## Be nice about error messages.
        if(n_of_objects(x) != n_of_objects(y))
            stop("Arguments 'x' and 'y' must have the same number of objects.")
        x <- cl_ensemble(x, y)
    }

    if(inherits(x, "cl_partition_ensemble"))
        .cl_join_partition(x)
    else if(inherits(x, "cl_dendrogram_ensemble"))
        .cl_join_dendrogram(x)
    else if(inherits(x, "cl_hierarchy_ensemble"))
        .cl_join_hierarchy(x)
    else
        stop("Cannot compute join of given clusterings.")
}

.cl_join_partition <-
function(x)
{
    x <- unique(x)
    if(length(x) == 1)
        return(cl_partition_by_class_ids(cl_class_ids(x[[1L]])))

    ## Canonicalize: ensure that class ids are always the integers from
    ## one to the number of classes.
    n <- sapply(x, n_of_classes)
    ids <- mapply(function(p, ncp) match(cl_class_ids(p),
                                         seq_len(ncp)),
                  x, n, SIMPLIFY = FALSE)
    ## Order according to the number of classes.
    ids <- ids[order(n)]

    ## And now incrementally build the join.
    jcids <- ids[[1L]]                  # Class ids of the current join.
    jnc <- length(unique(jcids))        # Number of classes of this.
    for(b in seq.int(from = 2, to = length(x))) {
        z <- table(jcids, ids[[b]])
        ## It is faster to work on the smaller partition, but this
        ## should be ensured by the reordering ...
        ## We need to "join all elements in the same class in at least
        ## one of the partitions".  In the matrix
        ##   C <- (tcrossprod(z) > 0)
        ## entry i,j is true/one iff z_{ik} z_{jk} > 0 for classes
        ## i and j in the current join (ids jcids) and some class k in
        ## the partition with ids[[b]], so that i and j must be joined.
        ## I.e., C indicates which classes need to be joined directly.
        ## We need to determine the transitive closure of this relation,
        ## which can be performed by repeating
        ##   C_{t+1} <- ((C_t %*% C) > 0)
        ## with C_1 = C until C_t does not change.
        C_new <- C_old <- C <- (tcrossprod(z) > 0)
        repeat {
            C_new <- (C_old %*% C) > 0
            if(all(C_new == C_old)) break
            C_old <- C_new
        }
        C <- C_new
        ## This should now have the connected components.
        ## Next, compute the map of the join class ids to the ids of
        ## these components. 
        cnt <- 0
        map <- remaining_ids <- seq_len(jnc)
        while(length(remaining_ids)) {
            cnt <- cnt + 1
            pos <- which(C[remaining_ids[1L], remaining_ids] > 0)
            map[remaining_ids[pos]] <- cnt
            remaining_ids <- remaining_ids[-pos]
        }
        ## And update the join:
        jcids <- map[jcids]
        jnc <- cnt
    }

    cl_partition_by_class_ids(jcids)    
}

.cl_join_dendrogram <-
function(x)
{
    ## Join of an ensemble of dendrograms.
    as.cl_dendrogram(do.call(pmax, lapply(x, cl_ultrametric)))
}

.cl_join_hierarchy <-
function(x)
{
    ## Join of an ensemble of n-trees.
    ## Only exists if the union of all classes of the n-trees is itself
    ## an n-tree (see Barthelemy et al).
    classes <- unique(unlist(lapply(x, cl_classes), recursive = FALSE))
    ## Now check if this is an n-tree.
    ## We must verify that for all classes A and B, their intersection
    ## is A, B, or empty.
    check <- function(A, B) {
        m_AB <- match(A, B)
        m_BA <- match(B, A)
        ((all(is.na(m_AB)) && all(is.na(m_BA)))
         || all(is.finite(m_AB))
         || all(is.finite(m_BA)))
    }
    for(i in seq_along(classes)) {
        A <- classes[[i]]
        for(j in seq_along(classes))
            if(!check(A, classes[[j]]))
                stop("Join of given n-trees does not exist.")
    }
    as.cl_hierarchy(.cl_ultrametric_from_classes(classes))
}
