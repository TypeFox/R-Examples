### * n_of_classes

## Get the number of classes in a (hard or soft) partition.

## <NOTE>
## We generally allow for classes to be empty, unlike the current
## version of kmeans().  Package cclust has a version of k-means which
## does not stop in case of empty classes.
## However, we only count NON-EMPTY classes here.
## </NOTE>

n_of_classes <-
function(x)
    UseMethod("n_of_classes")

## Default method.
n_of_classes.default <-
function(x)
    length(unique(cl_class_ids(x)))

## Package stats: kmeans() (R 2.1.0 or better).
n_of_classes.kmeans <- n_of_classes.default

## Package cluster: clara(), fanny(), and pam() give objects of the
## respective class inheriting from class "partition".
n_of_classes.fanny <-
function(x)
    sum(colSums(x$membership, na.rm = TRUE) > 0)
n_of_classes.partition <- n_of_classes.default

## Package cclust: cclust().
n_of_classes.cclust <- n_of_classes.default

## Package e1071: cmeans() gives objects of class "fclust".
n_of_classes.fclust <- n_of_classes.fanny
## Package e1071: cshell().
n_of_classes.cshell <- n_of_classes.fanny
## Package e1071: bclust().
n_of_classes.bclust <- n_of_classes.default

## Package mclust: Mclust().
n_of_classes.Mclust <- n_of_classes.default

## Package clue: Memberships.
n_of_classes.cl_membership <-
function(x)
    attr(x, "n_of_classes")
## Package clue: pclust().
n_of_classes.pclust <-
function(x)
{
    if(is.null(m <- x$membership)) 
        length(unique(cl_class_ids(x)))
    else
        sum(colSums(m, na.rm = TRUE) > 0)
}
## Package clue: (virtual) class "cl_partition".
n_of_classes.cl_partition <-
function(x)
    n_of_classes(.get_representation(x))

### * cl_class_ids

## Get ids of classes in a partition.
## <NOTE>
## Currently, all supported soft partitioning methods provide a softmax
## hard partitioning as well.
## </NOTE>

cl_class_ids <-
function(x)
    UseMethod("cl_class_ids")

## Default method.
cl_class_ids.default <-
function(x)
{
    stop("Cannot infer class ids from given object.")
}

## Package stats: kmeans() (R 2.1.0 or better).
cl_class_ids.kmeans <-
function(x)
    as.cl_class_ids(x$cluster)

## Package cluster: clara(), fanny(), and pam() give objects of the
## respective class inheriting from class "partition".
cl_class_ids.partition <-
function(x)
    as.cl_class_ids(x$clustering)

## Package RWeka: clusterers return objects inheriting from
## "Weka_clusterer".
cl_class_ids.Weka_clusterer <-
function(x)
    as.cl_class_ids(x$class_ids)

## Package cba: ccfkms().
cl_class_ids.ccfkms <-
function(x)
    as.cl_class_ids(as.vector(x$cl))
## Package cba: rockCluster() returns objects of class "rock".
cl_class_ids.rock <-
function(x)
    as.cl_class_ids(as.vector(x$cl))

## Package cclust: cclust().
cl_class_ids.cclust <- cl_class_ids.kmeans

## Package e1071: cmeans() gives objects of class "fclust".
cl_class_ids.fclust <- cl_class_ids.kmeans
## Package e1071: cshell().
cl_class_ids.cshell <- cl_class_ids.kmeans
## Package e1071: bclust().
cl_class_ids.bclust <- cl_class_ids.kmeans

## Package flexclust: kcca() returns objects of S4 class "kcca" which
## extends S4 class "flexclust".
## <NOTE>
## We used to be able to call flexclust::cluster(), but this now only
## has S4 methods for modeltools::clusters() S4 generic.  Let's call this
## one, and hope that flexclust has been loaded ...
## </NOTE>
cl_class_ids.kcca <-
function(x)
    as.cl_class_ids(modeltools::clusters(x))

## Package flexmix: class "flexmix".
## <NOTE>
## We used to be able to call flexmix::cluster(), but this now only has
## S4 methods for modeltools::clusters() S4 generic.  Let's call this
## one, and hope that flexmix has been loaded ...
## </NOTE>
cl_class_ids.flexmix <-
function(x)
    as.cl_class_ids(modeltools::clusters(x))

## Package kernlab: specc() and kkmeans() return objects of S4 class
## "specc".
cl_class_ids.specc <-
function(x)
{
    tmp <- unclass(x)
    as.cl_class_ids(.structure(as.vector(tmp), names = names(tmp)))
}

## Package mclust: Mclust().
cl_class_ids.Mclust <-
function(x)
    as.cl_class_ids(x$classification)

## Package relations: equivalence and preference relations.
cl_class_ids.relation <-
function(x)
    as.cl_class_ids(relations::relation_class_ids(x))

## Package clue: Class ids.
cl_class_ids.cl_class_ids <- identity
## Package clue: Memberships.
cl_class_ids.cl_membership <-
function(x)
    as.cl_class_ids(.structure(max.col(x), names = rownames(x)))
## (Cannot do cl_class_ids.cl_membership <- max.col for generic/method
## consistency.)
## Package clue: cl_pam().
cl_class_ids.cl_pam <- cl_class_ids.kmeans
## Package clue: cl_partition_by_class_ids().
cl_class_ids.cl_partition_by_class_ids <-
function(x)
    .get_representation(x)
## Package clue: kmedoids().
cl_class_ids.kmedoids <- cl_class_ids.kmeans
## Package clue: pclust().
cl_class_ids.pclust <- cl_class_ids.kmeans
## Package clue: (virtual) class "cl_partition".
cl_class_ids.cl_partition <-
function(x)
    cl_class_ids(.get_representation(x))

## Package movMF: class "movMF".
cl_class_ids.movMF <-
function(x)
    as.cl_class_ids(max.col(x$P))

### * as.cl_class_ids

as.cl_class_ids <-
function(x)
{
    ## For the time being, handle only "raw" class ids.
    ## Maybe add methods handling factors lateron (if necessary).
    ## <NOTE>
    ## This could also be used to canonicalize returned class ids
    ## according to the docs (vector of integers with the class ids),
    ## using someting like
    ##   match(ids, unique(ids))
    ## </NOTE>
    .structure(unclass(x), class = "cl_class_ids")
}

### * print.cl_class_ids

print.cl_class_ids <-
function(x, ...)
{
    writeLines("Class ids:")
    print(unclass(x), ...)
    invisible(x)
}

### * cl_class_labels

cl_class_labels <-
function(x)
    UseMethod("cl_class_labels")

### * is.cl_partition

## Determine whether an object is a (generalized) partition.
## Note that this includes both hard and soft partitions, and allows
## sums of memberships of objects to be less than one.

is.cl_partition <-
function(x)
    UseMethod("is.cl_partition")

## Default method.
is.cl_partition.default <- .false

## Package stats: kmeans() (R 2.1.0 or better).
is.cl_partition.kmeans <- .true

## Package cluster: clara(), fanny(), and pam() give objects of the
## respective class inheriting from class "partition".
is.cl_partition.partition <- .true

## Package RWeka: clusterers return objects inheriting from
## "Weka_clusterer".
## (Note that Cobweb internally uses a classification tree, but
## definitely does not expose this structure.)
is.cl_partition.Weka_clusterer <- .true

## Package cba: ccfkms().
is.cl_partition.ccfkms <- .true
## Package cba: rockCluster() returns objects of class "rock".
is.cl_partition.rock <- .true

## Package cclust: cclust().
is.cl_partition.cclust <- .true

## Package e1071: cmeans() gives objects of class "fclust".
is.cl_partition.fclust <- .true
## Package e1071: cshell().
is.cl_partition.cshell <- .true
## Package e1071: bclust().
is.cl_partition.bclust <- .true

## Package flexclust: kcca() returns objects of S4 class "kcca" which
## extends S4 class "flexclust".
is.cl_partition.kcca <- .true

## Package flexmix: class "flexmix".
is.cl_partition.flexmix <- .true

## Package kernlab: specc() and kkmeans() return objects of S4 class
## "specc".
is.cl_partition.specc <- .true

## Package mclust: Mclust().
is.cl_partition.Mclust <- .true

## Package clue: (virtual) class "cl_partition".
## Note that "raw" cl_membership objects are *not* partitions, as they
## are meant for numeric computations.
is.cl_partition.cl_partition <- .true
## Package clue: kmedoids().
is.cl_partition.kmedoids <- .true
## Package clue: pclust().
is.cl_partition.pclust <- .true

## Package movMF: class "movMF".
is.cl_partition.movMF <- .true

### * as.cl_partition

## Note that cl_partition conceptually is a virtual class, so there are
## no prototypes and no cl_partition() creator.

.cl_partition_classes <- "cl_partition"

as.cl_partition <-
function(x)
{
    if(is.cl_partition(x)) {
        if(!inherits(x, "cl_partition"))
            .make_container(x, .cl_partition_classes)
        else
            x
    }
    else
        cl_partition_by_memberships(as.cl_membership(x))
}

### * print.cl_partition

print.cl_partition <-
function(x, ...)
    .print_container(x, "cl_partition", ...)

### * print.cl_partition_by_class_ids

print.cl_partition_by_class_ids <-
function(x, ...)
{
    writeLines(gettextf("A hard partition of %d objects.",
                        n_of_objects(x)))
    print(cl_class_ids(x), ...)
    invisible(x)
}

### * print.cl_partition_by_memberships

print.cl_partition_by_memberships <-
function(x, ...)
{
    writeLines(gettextf("A partition of %d objects.",
                        n_of_objects(x)))
    print(cl_membership(x), ...)
    invisible(x)
}

### * Complex.cl_partition

Complex.cl_partition <-
function(z)
    stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                  .Generic, .Class),
         domain = NA)

### * Math.cl_partition

Math.cl_partition <-
function(x, ...)
    stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                  .Generic, .Class),
         domain = NA)

### * Ops.cl_partition

Ops.cl_partition <-
function(e1, e2)
{
    if(nargs() == 1L)
        stop(gettextf("Unary '%s' not defined for \"%s\" objects.",
                      .Generic, .Class),
             domain = NA)

    ## Only comparisons are supprorted.
    if(!(as.character(.Generic) %in% c("<", "<=", ">", ">=",
                                       "==", "!=")))
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class),
             domain = NA)

    ci1 <- cl_class_ids(e1)
    ci2 <- cl_class_ids(e2)
    if(length(ci1) != length(ci2))
        stop("Partitions must have the same number of objects.")
    z <- table(ci1, ci2) > 0
    switch(.Generic,
           "<=" = all(rowSums(z) == 1),
           "<"  = all(rowSums(z) == 1) && any(colSums(z) > 1),
           ">=" = all(colSums(z) == 1),
           ">"  = all(colSums(z) == 1) && any(rowSums(z) > 1),
           "==" = all(rowSums(z) == 1) && all(colSums(z) == 1),
           "!=" = any(rowSums(z) > 1) || any(colSums(z) > 1))

}

### * Summary.cl_partition

Summary.cl_partition <-
function(..., na.rm = FALSE)
{
    ok <- switch(.Generic, max = , min = , range = TRUE, FALSE)
    if(!ok)
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class),
             domain = NA)
    args <- list(...)
    switch(.Generic,
           "min" = cl_meet(cl_ensemble(list = args)),
           "max" = cl_join(cl_ensemble(list = args)),
           "range" = {
               cl_ensemble(min = cl_meet(cl_ensemble(list = args)),
                           max = cl_join(cl_ensemble(list = args)))
           })
}

### * cl_partition_by_class_ids

cl_partition_by_class_ids <-
function(x, labels = NULL)
{
    if(!is.atomic(x))
        stop("Class ids must be atomic.")
    if(is.null(names(x)))
        names(x) <- labels
    ## <FIXME>
    ## Perhaps give the raw class ids more structure?
    ## E.g, class "cl_class_ids"?
    ## Problem is that we used to say about extensibility that all there
    ## is to do for a hard partitioner is to add a cl_class_ids() method
    ## and two predicates, but *not* to have the former give a suitably
    ## classed object.  On the other hand, the recipe would need to be
    ## extended for soft partitioners, for which it would be necessary
    ## to provide a cl_membership() method which really returns an
    ## object of class cl_membership.  Note that we can do this using
    ## as.cl_membership(m), where m is the raw membership matrix.  So
    ## maybe we should ask for using as.cl_class_ids() to coerce raw
    ## class ids ...
    .make_container(as.cl_class_ids(x),
                    c("cl_partition_by_class_ids",
                      .cl_hard_partition_classes),
                    list(n_of_objects = length(x),
                         n_of_classes = length(unique(x))))
    ## </FIXME>
}

### * cl_partition_by_memberships

cl_partition_by_memberships <-
function(x, labels = NULL)
{
    if(!is.matrix(x)
       || any(x < 0, na.rm = TRUE)
       || any(x > 1, na.rm = TRUE))
        stop("Not a valid membership matrix.")
    ## Be nice.
    x <- x / rowSums(x, na.rm = TRUE)
    ## (Note that this does not imply all(rowSums(x) == 1).  If we
    ## wanted to test for this, something like
    ##    .is_stochastic_matrix <- function(x)
    ##       identical(all.equal(rowSums(x), rep(1, nrow(x))), TRUE))
    ## should do.)
    if(is.null(rownames(x)))
        rownames(x) <- labels
    .make_container(as.cl_membership(x),
                    c("cl_partition_by_memberships",
                      .cl_partition_classes),
                    list(n_of_objects = nrow(x)))
}

### * is.cl_hard_partition

## Determine whether an object is a hard partition.

is.cl_hard_partition <-
function(x)
    UseMethod("is.cl_hard_partition")

## Default method.
is.cl_hard_partition.default <- .false

## Package stats: kmeans() (R 2.1.0 or better).
is.cl_hard_partition.kmeans <- .true

## Package cluster: clara(), fanny(), and pam() give objects of the
## respective class inheriting from class "partition".
## <NOTE>
## Of course, fuzzy clustering can also give a hard partition ...
is.cl_hard_partition.fanny <-
function(x)
{
    all(rowSums(cl_membership(x) == 1, na.rm = TRUE) > 0)
}
## </NOTE>
is.cl_hard_partition.partition <- .true

## Package RWeka: clusterers return objects inheriting from
## "Weka_clusterer".
is.cl_hard_partition.Weka_clusterer <- .true

## Package cba: ccfkms().
is.cl_hard_partition.ccfkms <- .true
## Package cba: rockCluster() returns objects of class "rock".
is.cl_hard_partition.rock <- .true

## Package cclust: cclust().
is.cl_hard_partition.cclust <- .true

## Package e1071: cmeans() gives objects of class "fclust".
is.cl_hard_partition.fclust <- is.cl_hard_partition.fanny
## Package e1071: cshell().
is.cl_hard_partition.cshell <- is.cl_hard_partition.fanny
## Package e1071: bclust().
is.cl_hard_partition.bclust <- .true

## Package flexclust: kcca() returns objects of S4 class "kcca" which
## extends S4 class "flexclust".
is.cl_hard_partition.kcca <- .true

## Package flexmix: class "flexmix".
is.cl_hard_partition.flexmix <- is.cl_hard_partition.fanny

## Package kernlab: specc() and kkmeans() return objects of S4 class
## "specc".
is.cl_hard_partition.specc <- .true

## Package mclust: Mclust().
is.cl_hard_partition.Mclust <- is.cl_hard_partition.fanny

## Package clue: (virtual) class "cl_hard_partition".
is.cl_hard_partition.cl_hard_partition <- .true
## Package clue: (virtual) class "cl_partition".
## Note that "raw" cl_membership objects are *not* partitions, as they
## are meant for numeric computations.
## Rather than providing is.cl_hard_partition.cl_membership() we thus
## prefer explicit handling of cl_partition objects with a cl_membership
## representation.
is.cl_hard_partition.cl_partition <-
function(x)
{
    ## If the object has a cl_membership representation ...
    y <- .get_representation(x)
    if(inherits(y, "cl_membership"))
        attr(y, "is_cl_hard_partition")
    ## Other representations, e.g. for "definitely" hard partitions via
    ## vectors of class ids or class labels, or a list of classes, may
    ## be added in future versions.
    ## In any case, this must be kept in sync with what is handled by
    ## as.cl_partition() [which currently runs as.cl_membership() in
    ## case is.cl_partition() gives false].
    else
        is.cl_hard_partition(y)
}
## Package clue: kmedoids().
is.cl_hard_partition.kmedoids <- .true
## Package clue: pclust().
is.cl_hard_partition.pclust <- is.cl_hard_partition.fanny

## Package movMF: class "movMF".
is.cl_hard_partition.movMF <- is.cl_hard_partition.fanny

### * as.cl_hard_partition

.cl_hard_partition_classes <- c("cl_hard_partition", "cl_partition")

as.cl_hard_partition <-
function(x)
{
    if(is.cl_hard_partition(x)) {
        if(!inherits(x, "cl_partition"))
            .make_container(x, .cl_hard_partition_classes)
        else
            x
    }
    else if(is.cl_partition(x)) {
        ## A soft cl_partition ...
        ids <- cl_class_ids(x)
        cl_partition_by_class_ids(ids, names(ids))
    }
    else if(is.matrix(x)) {
        ## A matrix of raw memberships, hopefully ...
        cl_partition_by_class_ids(max.col(x), rownames(x))
    }
    else if(is.atomic(x)) {
        ## A vector of raw class ids, hopefully ...
        cl_partition_by_class_ids(x, names(x))
    }
    else
        stop("Cannot coerce to 'cl_hard_partition'.")
}

### * is.cl_soft_partition

## Determine whether an object is a soft partition.

is.cl_soft_partition <-
function(x)
    is.cl_partition(x) && ! is.cl_hard_partition(x)

### * .maybe_is_proper_soft_partition

## Determine whether an object might be a proper soft partition (in the
## sense that it is a cl_partition but not a cl_hard_partition).
## This is mostly useful when computing fuzziness measures.

.maybe_is_proper_soft_partition <-
function(x)
    UseMethod(".maybe_is_proper_soft_partition")
.maybe_is_proper_soft_partition.default <- .false
.maybe_is_proper_soft_partition.fanny <- .true
.maybe_is_proper_soft_partition.fclust <- .true
.maybe_is_proper_soft_partition.cshell <- .true
.maybe_is_proper_soft_partition.flexmix <- .true
.maybe_is_proper_soft_partition.Mclust <- .true
## See above for why we prefer not to have
## .maybe_is_proper_soft_partition.cl_membership().
## (Although this is an internal generic really only used for making
## cl_fuzziness() computations more efficient, so we could be more
## generous here [perhaps using a slightly different name such as
## .maybe_represents_a_proper_soft_partition()].
.maybe_is_proper_soft_partition.cl_partition <-
function(x)
{
    y <- .get_representation(x)
    if(inherits(y, "cl_membership"))
        !attr(y, "is_cl_hard_partition")
    else
        .maybe_is_proper_soft_partition(y)
}
.maybe_is_proper_soft_partition.pclust <-
function(x)
    x$m > 1


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
