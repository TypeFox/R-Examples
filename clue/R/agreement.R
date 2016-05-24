### * cl_agreement

cl_agreement <-
function(x, y = NULL, method = "euclidean", ...)
{
    ## <NOTE>
    ## This code is repeated from cl_dissimilarity(), mutatis mutandis.
    ## Not really a big surprise ...
    ## </NOTE>
    
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
        if(!inherits(method, "cl_agreement_method")) {
            ## Get the method definition and description from the
            ## registry. 
            type <- ifelse(is_partition_ensemble,
                           "partition", "hierarchy")
            method <- get_cl_agreement_method(method, type)
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
        ## Build a cross-proximity object of cross-agreements.
        d <- matrix(0, length(x), length(y))
        for(j in seq_along(y))
            d[, j] <- sapply(x, method, y[[j]], ...)
        dimnames(d) <- list(names(x), names(y))
        return(cl_cross_proximity(d, method_name,
                                  class = "cl_cross_agreement"))
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

    ## <NOTE>
    ## We assume that self-agreements are always one ...
    ## </NOTE>
    cl_proximity(unlist(d),
                 method_name,
                 labels = names(x),
                 self = rep.int(1, length(x)),
                 size = n, class = "cl_agreement")
}

### ** .cl_agreement_partition_euclidean

.cl_agreement_partition_euclidean <-
function(x, y)
{
    ## <NOTE>
    ## Upper bound for maximal dissimilarity, maybe improve eventually.
    d_max <- sqrt(2 * n_of_objects(x))
    ## </NOTE>
    1 - .cl_dissimilarity_partition_euclidean(x, y) / d_max
}

### ** .cl_agreement_partition_manhattan

.cl_agreement_partition_manhattan <-
function(x, y)
{
    ## <NOTE>
    ## Upper bound for maximal dissimilarity, maybe improve eventually.
    d_max <- 2 * n_of_objects(x)
    ## </NOTE>
    1 - .cl_dissimilarity_partition_manhattan(x, y) / d_max
}
    
### ** .cl_agreement_partition_Rand

.cl_agreement_partition_Rand <-
function(x, y)
{
    n <- n_of_objects(x)
    
    ## Handle soft partitions using the corresponding hard ones.
    ## (At least, for the time being.)
    x <- table(cl_class_ids(x), cl_class_ids(y))

    ## <NOTE>
    ## The number A of concordant pairs is given by
    ##   A = choose(n,2) + \sum_{i,j} x_{ij}^2
    ##       - (1/2) * (\sum_i x_{i.}^2 + \sum_j x_{.j}^2)
    ##     = choose(n,2) + 2 \sum_{i,j} choose(x_{ij},2)
    ##       - (\sum_i choose(x_{i.},2) + \sum_j choose(x_{.j},2)
    ## with the first version certainly much faster to compute.
    ## </NOTE>
    1 + (sum(x^2) -
         (sum(rowSums(x)^2) + sum(colSums(x)^2)) / 2) / choose(n, 2)
}

### ** .cl_agreement_partition_cRand

.cl_agreement_partition_cRand <-
function(x, y)
{
    if(!is.cl_hard_partition(x) || !is.cl_hard_partition(y))
        stop("Can only handle hard partitions.")

    n <- n_of_objects(x)    
    x <- table(cl_class_ids(x), cl_class_ids(y))

    ## <NOTE>
    ## The basic formula is
    ##   (Sxy - E) / ((Sx. + S.y) / 2 - E)
    ## where
    ##   Sxy = \sum_{i,j} choose(x_{ij}, 2)
    ##   Sx. = \sum_i     choose(x_{i.}, 2)
    ##   S.y = \sum_j     choose(x_{.j}, 2)
    ## and
    ##   E = Sx. * S.y / choose(n, 2)
    ## We replace the bincoefs by the corresponding sums of squares,
    ## getting
    ##   (Txy - F) / ((Tx. + T.y) / 2 - F)
    ## where
    ##   Txy = \sum_{i,j} x_{ij}^2 - n
    ##   Tx. = \sum_i     x_{i.}^2 - n
    ##   T.y = \sum_j     x_{.j}^2 - n
    ## and
    ##   F = Tx. * T.y / (n^2 - n)
    ## </NOTE>
    Txy <- sum(x ^ 2) - n
    Tx. <- sum(rowSums(x) ^ 2) - n
    T.y <- sum(colSums(x) ^ 2) - n
    F <- Tx. * T.y / (n ^ 2 - n)
    (Txy - F) / ((Tx. + T.y) / 2 - F)
}


### ** .cl_agreement_partition_NMI

.cl_agreement_partition_NMI <-
function(x, y)
{
    if(!is.cl_hard_partition(x) || !is.cl_hard_partition(y))
        stop("Can only handle hard partitions.")
    x <- table(cl_class_ids(x), cl_class_ids(y))
    x <- x / sum(x)
    
    m_x <- rowSums(x)
    m_y <- colSums(x)
    y <- outer(m_x, m_y)

    i <- which((x > 0) & (y > 0))
    out <- sum(x[i] * log(x[i] / y[i]))
    e_x <- sum(m_x * log(ifelse(m_x > 0, m_x, 1)))
    e_y <- sum(m_y * log(ifelse(m_y > 0, m_y, 1)))

    out / sqrt(e_x * e_y)
}

### ** .cl_agreement_partition_KP

.cl_agreement_partition_KP <-
function(x, y)
{
    ## Agreement measure due to Katz & Powell (1953, Psychometrika), see
    ## also Messatfa (1992, Journal of Classification).

    n <- n_of_objects(x)    

    ## Handle soft partitions using the corresponding hard ones.
    ## (At least, for the time being.)
    x <- table(cl_class_ids(x), cl_class_ids(y))

    A_xy <- sum(x ^ 2)
    A_x. <- sum(rowSums(x) ^ 2)
    A_.y <- sum(colSums(x) ^ 2)

    (n^2 * A_xy - A_x. * A_.y) /
        sqrt(A_x. * (n^2 - A_x.) * A_.y * (n^2 - A_.y))
}

### ** .cl_agreement_partition_angle

.cl_agreement_partition_angle <-
function(x, y)
{
    ## Maximal angle between the matched memberships.

    k <- max(n_of_classes(x), n_of_classes(y))
    M_x <- cl_membership(x, k)
    M_y <- cl_membership(y, k)
    ## Match classes from conforming memberships.
    ind <- solve_LSAP(crossprod(M_x, M_y), maximum = TRUE)
    sum(M_x * M_y[, ind]) / sqrt(sum(M_x ^ 2) * sum(M_y ^ 2))
}

### ** .cl_agreement_partition_diag

.cl_agreement_partition_diag <-
function(x, y)
{
    ## Maximal co-classification rate.

    k <- max(n_of_classes(x), n_of_classes(y))
    M_x <- cl_membership(x, k)
    M_y <- cl_membership(y, k)
    ## Match classes from conforming memberships.
    ind <- solve_LSAP(crossprod(M_x, M_y), maximum = TRUE)
    sum(M_x * M_y[, ind]) / n_of_objects(x)
}

### ** .cl_agreement_partition_FM

.cl_agreement_partition_FM <-
function(x, y)
{
    ## Fowlkes-Mallows index.

    n <- n_of_objects(x)    

    ## Handle soft partitions using the corresponding hard ones.
    ## (At least, for the time being.)
    x <- table(cl_class_ids(x), cl_class_ids(y))

    (sum(x ^ 2) - n) /
        sqrt((sum(rowSums(x) ^ 2) - n) * (sum(colSums(x) ^ 2) - n))
    
}

### ** .cl_agreement_partition_Jaccard

.cl_agreement_partition_Jaccard <-
function(x, y)
{
    ## Jaccard index.

    n <- n_of_objects(x)    

    ## Handle soft partitions using the corresponding hard ones.
    ## (At least, for the time being.)
    x <- table(cl_class_ids(x), cl_class_ids(y))

    Z <- sum(x ^ 2)

    (Z - n) / (sum(rowSums(x) ^ 2) + sum(colSums(x) ^ 2) - n - Z)
    
}

### ** .cl_agreement_partition_purity

.cl_agreement_partition_purity <-
function(x, y)
{
    ## Purity of classes of x with respect to those of y: relative
    ## fraction of "optimally matched and collapsed" joint class
    ## frequencies, i.e., \sum_i \max_j c_{ij} / n.

    n <- n_of_objects(x)    

    ## Handle soft partitions using the corresponding hard ones.
    ## (At least, for the time being.)
    x <- table(cl_class_ids(x), cl_class_ids(y))

    sum(apply(x, 1L, max)) / n
}

.cl_agreement_partition_PS <-
function(x, y)
{
    ## Prediction Strength as used in Tibshirani and Walter (2005),
    ## "Cluster Validation by Prediction Strength", JCGS.

    ## See Eqn 2.1 in the reference: this is
    ##   min_l rate of different objects in the same class in partition
    ##         A and in class l in partition B,
    ## where the min is taken over all classes l of partition B.

    x <- table(cl_class_ids(x), cl_class_ids(y))
    s <- rowSums(x)

    min((rowSums(x ^ 2) - s) / (s * (s - 1)), na.rm = TRUE)
}   

## Some computations useful for interpreting some of the above.
##
## Consider two hard partitions A and B and write
##   a_{ik} ... indicator of object i in class k for partition A
##   b_{il} ... indicator of object i in class l for partition B
## (so that the a_{ik} and b_{il} are of course the membership matrices
## of the partitions).
##
## Then obviously
##   \sum_i a_{ik} b_{il} = m_{kl}
## is the number of objects in class k for A and in class l for B, and
##   \sum_i a_{ik} = m_{k.} = # objects in class k for A
##   \sum_i b_{il} = m_{.l} = # objects in class l for B
##
## Number of pairs of objects in the same classes for both A and B:
##   \sum_{i, j, k, l} a_{ik} a_{jk} b_{il} b_{jl}
##     = \sum_{k, l} \sum_i a_{ik} b_{il} \sum_j a_{jk} b_{jl}
##     = \sum_{k, l} m_{kl} ^ 2
## This includes the n pairs with identical objects, hence:
## Number of distinct pairs of objects in the same classes for both A
## and B:
##   (\sum_{k, l} m_{kl} ^ 2 - n) / 2
##
## Number of pairs of objects in the same class for A:
##   \sum_{i, j, k} a_{ik} a_{jk}
##     = \sum_k \sum_i a_{ik} \sum_j a_{jk}
##     = \sum_k m_{k.} ^ 2
## Again, this includes the n pairs with identical objects, hence:
## Number of distinct pairs of objects in the same class for A:
##   (\sum_k m_{k.} ^ 2 - n) / 2
##
## Similarly, \sum_l m_{.l} ^ 2 corresponds to the number of pairs of
## objects in the same class for B.
##
## Finally, to get the number of pairs of objects in different classes
## for both A and B, we note that this is the total number of pairs,
## minus the sum of the numbers of those in the same class for A and for
## B, respectively, plus the number of pairs in the same class for both
## A and B.
##
## This makes e.g. the interpretation of some of the Fowlkes-Mallows or
## Rand agreement indices rather straightforward.

### ** .cl_agreement_hierarchy_euclidean

.cl_agreement_hierarchy_euclidean <-
function(x, y)
    1 / (1 + .cl_dissimilarity_hierarchy_euclidean(x, y))

### ** .cl_agreement_hierarchy_manhattan

.cl_agreement_hierarchy_manhattan <-
function(x, y)
    1 / (1 + .cl_dissimilarity_hierarchy_manhattan(x, y))

### ** .cl_agreement_hierarchy_cophenetic

.cl_agreement_hierarchy_cophenetic <-
function(x, y)
{
    ## Cophenetic correlation.
    if(!.has_object_dissimilarities(x) ||
       !.has_object_dissimilarities(y))
        return(NA)
    cor(cl_object_dissimilarities(x), cl_object_dissimilarities(y))
}

### ** .cl_agreement_hierarchy_angle

.cl_agreement_hierarchy_angle <-
function(x, y)
{
    ## Angle between ultrametrics.
    if(!.has_object_dissimilarities(x) ||
       !.has_object_dissimilarities(y))
        return(NA)
    u_x <- cl_object_dissimilarities(x)
    u_y <- cl_object_dissimilarities(y)
    sum(u_x * u_y) / sqrt(sum(u_x ^ 2) * sum(u_y ^ 2))
}

### ** .cl_agreement_hierarchy_gamma

.cl_agreement_hierarchy_gamma <-
function(x, y)
    1 - .cl_dissimilarity_hierarchy_gamma(x, y)
    
### * [.cl_agreement

"[.cl_agreement" <-
function(x, i, j)
{
    y <- NextMethod("[")
    if(!inherits(y, "cl_agreement")) {
        description <- attr(x, "description")
        return(cl_cross_proximity(y, description = description,
                                  class = "cl_cross_agreement"))
    }
    y
}


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
