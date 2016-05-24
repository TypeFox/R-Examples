### * cl_medoid

cl_medoid <-
function(x, method = "euclidean")
{
    ## <NOTE>
    ## In principle we can get the same using pam(k = 1)$medoids.
    ## </NOTE>

    clusterings <- as.cl_ensemble(x)
    if(!length(clusterings))
        stop("Cannot compute medoid of empty ensemble.")
    dissimilarities <-
        as.matrix(cl_dissimilarity(clusterings, method = method))
    clusterings[[which.min(rowSums(dissimilarities))]]
}

### * cl_pam

cl_pam <-
function(x, k, method = "euclidean", solver = c("pam", "kmedoids"))
{

    clusterings <- as.cl_ensemble(x)
    if(!length(clusterings))
        stop("Cannot compute medoid partition of empty ensemble.")

    ## Actually, we should have at least k distinct elements in the
    ## ensemble ...

    make_cl_pam <- function(class_ids, medoid_ids, medoids, criterion,
                            description)
        .structure(list(cluster = class_ids,
                        medoid_ids = medoid_ids,
                        prototypes = medoids,
                        criterion = criterion,
                        description = description),
                   class = "cl_pam")

    if(k == 1L) {
        ## Simplify matters if a global medoid is sought.
        dissimilarities <-
            cl_dissimilarity(clusterings, method = method)
        description <- attr(dissimilarities, "description")
        dissimilarities <- as.matrix(dissimilarities)
        row_sums <- rowSums(dissimilarities)
        medoid_id <- which.min(row_sums)
        criterion <- row_sums[medoid_id]
        return(make_cl_pam(as.cl_class_ids(seq_along(clusterings)),
                           medoid_id,
                           clusterings[medoid_id],
                           criterion,
                           description))
    }

    solver <- match.arg(solver)

    ## Argh.  We really want to run k-medoids for the unique elements of
    ## the ensemble, but pam() only works for symmetric dissimilarties.
    ## As computing cluster dissimilarities is typically expensive, use
    ## the unique elements for doing so in any case.
    
    values <- unique(clusterings)
    ## Positions of ensemble members in the unique values.
    positions <- match(clusterings, values)

    ## Dissimilarities between unique values.
    dissimilarities <- cl_dissimilarity(values, method = method)
    description <- attr(dissimilarities, "description")
    dissimilarities <- as.matrix(dissimilarities)

    ## For pam(), we need the dissimilarities for all objects.
    if(solver == "pam") {
        dissimilarities <- dissimilarities[positions, positions]
        party <- cluster::pam(as.dist(dissimilarities), k)
        class_ids <- cl_class_ids(party)
        medoid_ids <- cl_medoid_ids(party)
        medoids <- clusterings[medoid_ids]
        criterion <- sum(dissimilarities[cbind(seq_along(class_ids),
                                               medoid_ids[class_ids])])
    }
    else {
        ## Counts of unique values.
        counts <- tabulate(positions)

        ## Weigh according to the counts.  Should be straightforward to
        ## add "case weights" as well ...
        dissimilarities <- counts * dissimilarities

        ## Now partition.
        party <- kmedoids(dissimilarities, k)

        ## And build the solution from this ...
        criterion <- party$criterion        
        ## First, things for the unique values.
        medoid_ids <- cl_medoid_ids(party)
        medoids <- values[medoid_ids]
        class_ids <- cl_class_ids(party)
        ## Second, things for all objects.
        class_ids <- class_ids[positions]
        medoid_ids <- match(medoid_ids, positions)
    }

    make_cl_pam(class_ids, medoid_ids, medoids, criterion, description)
}

print.cl_pam <-
function(x, ...)
{
    class_ids <- cl_class_ids(x)
    fmt <- "A k-medoid partition of a cluster ensemble with %d elements into %d
classes (dissimilarity measure: %s)."
    writeLines(c(strwrap(gettextf(fmt,
                                  n_of_objects(x),
                                  n_of_classes(x),
                                  x$description))))
    writeLines(gettext("Class ids:"))
    print(class_ids, ...)
    writeLines(gettext("Criterion:"))
    print(x$criterion, ...)
    invisible(x)
}

### * cl_medoid_ids

## Little helper, internal for the time being ...
cl_medoid_ids <- function(x) UseMethod("cl_medoid_ids")
cl_medoid_ids.cl_pam <- function(x) x$medoid_ids
cl_medoid_ids.kmedoids <- function(x) x$medoid_ids
cl_medoid_ids.clara <- function(x) x$i.med
cl_medoid_ids.pam <- function(x) x$id.med

### * kmedoids

kmedoids <-
function(x, k)
{
    ## <FIXME>
    ## For the time being, 'x' is assumed a dissimilarity object or a
    ## matrix of dissimilarities.
    ## Let's worry about the interface later.
    ## </FIXME>

    x <- as.matrix(x)

    n <- nrow(x)

    ## Use the formulation in Gordon & Vichi (1998), Journal of
    ## Classification, [P4'], page 279, with variables c(vec(X), z), but
    ## with rows and cols interchanged (such that x_{ij} is one iff o_i
    ## has medoid o_j, and z_j is one iff o_j is a medoid).

    make_constraint_mat <- function(n) {
        nsq <- n * n
        rbind(cbind(kronecker(rbind(rep.int(1, n)), diag(1, n)),
                    matrix(0, n, n)),
              cbind(diag(1, nsq),
                    kronecker(diag(1, n), rep.int(-1, n))),
              c(double(nsq), rep.int(1, n)),
              cbind(matrix(0, n, nsq), diag(1, n)))
    }

    make_constraint_dir <- function(n)
        rep.int(c("=", "<=", "=", "<="), c(n, n * n, 1, n))

    make_constraint_rhs <- function(n, k)
        rep.int(c(1, 0, k, 1), c(n, n * n, 1, n))

    ## <NOTE>
    ## We could try a relaxation without integrality constraints first,
    ## which seems to "typically work" (and should be faster).  To test
    ## for integrality, use something like
    ##   if(identical(all.equal(y$solution, round(y$solution)), TRUE))
    ## </NOTE>

    y <- lpSolve::lp("min",
                     c(c(x), double(n)),
                     make_constraint_mat(n),
                     make_constraint_dir(n),
                     make_constraint_rhs(n, k),
                     int.vec = seq_len(n * (n + 1)))

    ## Now get the class ids and medoids.

    ind <- which(matrix(y$solution[seq_len(n * n)], n) > 0,
                 arr.ind = TRUE)
    medoid_ids <- unique(ind[, 2L])
    class_ids <- seq_len(n)
    class_ids[ind[, 1L]] <- match(ind[, 2L], medoid_ids)

    .structure(list(cluster = class_ids,
                    medoid_ids = medoid_ids,
                    criterion = y$objval),
               class = "kmedoids")
}

print.kmedoids <-
function(x, ...)
{
    fmt <- "A k-medoids clustering of %d objects into %d clusters."
    writeLines(gettextf(fmt, n_of_objects(x), n_of_classes(x)))
    writeLines(gettext("Medoid ids:"))
    print(cl_medoid_ids(x), ...)
    writeLines(gettext("Class ids:"))
    print(unclass(cl_class_ids(x)), ...)
    writeLines(gettext("Criterion:"))
    print(x$criterion, ...)
    invisible(x)
}




