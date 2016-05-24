### * cl_membership

## Get the class membership matrix from a partition.

## <NOTE>
## We could use sparse matrices for the memberships of hard partitions.
## Not sure if this is really that important, though, as we typically
## use memberships in a context where dense matrices (memberships of
## soft partitions) occur.
## </NOTE>

## <NOTE>
## Currently, the number of classes to be used for the memberships must
## not be less than the number of classes in the partition.  We might
## eventually change this so that "optimal" collapsing of classes is
## performed (but note that optimality needs to be relative to some
## dissimilarity measure) ...
## However, from the discussion of the second method in Gordon and Vichi
## (2001) we note that whereas optimal assignment is "simple", optimal
## collapsing (equivalent to partitioning into an arbitrary number of
## partitions) is of course very hard.
## </NOTE>

cl_membership <-
function(x, k = n_of_classes(x))
{
    if(k < n_of_classes(x))
        stop("k cannot be less than the number of classes in x.")
    UseMethod("cl_membership")
}

## Default method.
cl_membership.default <-
function(x, k = n_of_classes(x))
    .cl_membership_from_class_ids(cl_class_ids(x), k)

## Package stats: kmeans() (R 2.1.0 or better).
cl_membership.kmeans <- cl_membership.default

## Package cluster: clara(), fanny(), and pam() give objects of the
## respective class inheriting from class "partition".
cl_membership.fanny <-
function(x, k = n_of_classes(x))
    .cl_membership_from_memberships(x$membership, k)
cl_membership.partition <- cl_membership.default

## Package cclust: cclust().
cl_membership.cclust <- cl_membership.default

## Package e1071: cmeans() gives objects of class "fclust".
cl_membership.fclust <- cl_membership.fanny
## Package e1071: cshell().
cl_membership.cshell <- cl_membership.fanny
## Package e1071: bclust().
cl_membership.bclust <- cl_membership.default

## Package flexmix: class "flexmix".
## <NOTE>
## We used to be able to call flexmix::posterior(), but this now only
## has S4 methods for modeltools::posterior() S4 generic.  Let's call
## this one, and hope that flexmix has been loaded ...
## </NOTE>
cl_membership.flexmix <-
function(x, k = n_of_classes(x))
    .cl_membership_from_memberships(modeltools::posterior(x), k)

## Package mclust: Mclust().
cl_membership.Mclust <-
function(x, k = n_of_classes(x))
    .cl_membership_from_memberships(x$z, k)

## Package clue: Memberships.
cl_membership.cl_membership <-
function(x, k = n_of_classes(x))
    .cl_membership_from_memberships(x, k)
## (Note: we cannot simply return x in case k equals n_of_classes(x),
## because ncol(x) might be different.)

## Package clue: pclust().
cl_membership.pclust <-
function(x, k = n_of_classes(x))
{
    ## We should really have a suitable "sparse matrix" class for
    ## representing the memberships of hard partitions.  In case we
    ## decide not to fill the membership "slot" for such:
    if(is.null(m <- x$membership))
        .cl_membership_from_class_ids(x$cluster, k)
    else
        .cl_membership_from_memberships(m, k)
}

## Package clue: (virtual) class "cl_partition".
cl_membership.cl_partition <-
function(x, k = n_of_classes(x))
    cl_membership(.get_representation(x), k)

## Package movMF: class "movMF".
cl_membership.movMF <-
function(x, k = n_of_classes(x))
    .cl_membership_from_memberships(x$P, k)

### * .make_cl_membership

## A low-level common creator.

.make_cl_membership <-
function(x, n_of_classes, is_cl_hard_partition, meta = NULL)
{
    attr(x, "n_of_classes") <- n_of_classes
    attr(x, "is_cl_hard_partition") <- is_cl_hard_partition
    attr(x, "meta") <- meta
    class(x) <- "cl_membership"
    x
}

### * .cl_membership_from_class_ids

.cl_membership_from_class_ids <-
function(x, k = NULL, meta = NULL)
{
    x <- factor(x)
    n_of_objects <- length(x)
    n_of_classes <- nlevels(x)
    if(is.null(k))
        k <- n_of_classes
    else if(k < n_of_classes)
        stop("k cannot be less than the number of classes in x.")
    ## <TODO>
    ## Should really use a sparse encoding of this ...
    M <- matrix(0, n_of_objects, k)
    ## (Could also use .one_entry_per_column(M, as.numeric(x)) <- 1 for
    ## the time being.)
    M[cbind(seq_len(n_of_objects), as.numeric(x))] <- 1
    ## But note that we also need to handle NAs ...
    M[is.na(x), ] <- NA
    ## </TODO>
    if(nlevels(x) == k)
        colnames(M) <- levels(x)
    if(!is.null(nm <- names(x)))
        rownames(M) <- nm
    .make_cl_membership(M, n_of_classes, TRUE, meta)
}

### * .cl_membership_from_memberships

.cl_membership_from_memberships <-
function(x, k = NULL, meta = NULL)
{
    ## <NOTE>
    ## Dropping and re-filling of ## zero columns in case k is given may
    ## seem unnecessary, but really canonicalizes by moving zero columns
    ## last ...
    ## </NOTE>
    
    x <- x[ , colSums(x, na.rm = TRUE) > 0, drop = FALSE]
    n_of_classes <- ncol(x)
    if(!is.null(k)) {
        if(k < n_of_classes)
            stop("k cannot be less than the number of classes in x.")
        if(k > n_of_classes) {
            ## Fill up with zero columns.
            x <- cbind(x, matrix(0, nrow(x), k - n_of_classes))
            ## Handle NAs if necessary.
            x[apply(is.na(x), 1, any), ] <- NA
        }
    }
    .make_cl_membership(x, n_of_classes,
                        all(rowSums(x == 1, na.rm = TRUE) > 0),
                        meta)
}

### * as.cl_membership

as.cl_membership <-
function(x)
    UseMethod("as.cl_membership")
as.cl_membership.default <-
function(x)
{
    if(inherits(x, "cl_membership"))
        x
    else if(is.atomic(x))
        .cl_membership_from_class_ids(x)
    else
        cl_membership(x)
}
as.cl_membership.matrix <-
function(x)
    .cl_membership_from_memberships(x)

### * .memberships_from_cross_dissimilarities

.memberships_from_cross_dissimilarities <-
function(d, power = 2)
{
    ## For a given matrix of cross-dissimilarities [d_{bj}], return a
    ## matrix [u_{bj}] such that \sum_{b,j} u_{bj}^p d_{bj}^q => min!
    ## under the constraint that u is a stochastic matrix.
    ## If only one power is given, it is taken as p, with q as 1.
    ## <NOTE>
    ## This returns a plain matrix of membership values and not a
    ## cl_membership object (so that it does not deal with possibly
    ## dropping or re-introducing unused classes).
    ## </NOTE>
    exponent <- if(length(power) == 1L)
        1 / (1 - power)
    else
        power[2L] / (1 - power[1L])
    u <- matrix(0, nrow(d), ncol(d))
    zero_incidences <- !(d > 0)
    n_of_zeroes <- rowSums(zero_incidences)
    if(any(ind <- (n_of_zeroes > 0)))
        u[ind, ] <-
            zero_incidences[ind, , drop = FALSE] / n_of_zeroes[ind]
    if(any(!ind)) {
        ## Compute d_{bj}^e / \sum_k d_{bk}^e without overflow from very
        ## small d_{bj} values.
        d <- exponent * log(d[!ind, , drop = FALSE])
        d <- exp(d - d[cbind(seq_len(nrow(d)), max.col(d))])
        u[!ind, ] <- d / rowSums(d)
    }
    u
}

### * print.cl_membership

print.cl_membership <-
function(x, ...)
{
    writeLines("Memberships:")
    print(matrix(as.vector(x), nrow = nrow(x), dimnames = dimnames(x)),
          ...)
    invisible(x)
}

### .has_object_memberships

## Be nice to users when computing proximities: all measures for
## "partitions" we currently consider really only assume that we can
## compute memberships and/or class ids.

## Note that the cl_membership() default method works for cl_class_ids.

.has_object_memberships <-
function(x)
    (is.cl_partition(x)
     || inherits(x, "cl_membership")
     || inherits(x, "cl_class_ids"))

### * .stochastify

.stochastify <-
function(x)
{
    ## Try to ensure that a stochastic matrix is returned.
    x <- pmax(x, 0)
    x / rowSums(x)
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
