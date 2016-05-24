## A slightly polymorphic function, similar to cluster::silhouette() and
## its methods.

cl_validity <-
function(x, ...)
    UseMethod("cl_validity")

cl_validity.default <-
function(x, d, ...)
{
    ## Note that providing methods for classes "cl_partition" and
    ## "cl_hierarchy" is not good enough ...
    out <- list()
    if(.has_object_memberships(x)) {
        v <- .cl_validity_partition_d_a_f(cl_membership(x),
                                          as.matrix(d))
        out <- list("Dissimilarity accounted for" = v)
    }
    else if(.has_object_dissimilarities(x)) {
        x <- cl_object_dissimilarities(x)
        d <- as.dist(d)
        out <- list("Variance accounted for" =
                    .cl_validity_hierarchy_variance_a_f(x, d),
                    "Deviance accounted for" =
                    .cl_validity_hierarchy_deviance_a_f(x, d))
        ## Consider adding e.g. the Agglomerative Coefficient or
        ## Divisive Coeffcient for more than cluster::agnes() and
        ## cluster::diana(), respectively.
    }
    class(out) <- "cl_validity"
    out
}

## Package cluster: agnes().
cl_validity.agnes <-
function(x, ...)
{
    out <- list("Agglomerative coefficient" = x$ac)
    ## According to the docs, agnes objects always have a diss
    ## component, but let's be defensive ...
    if(!is.null(d <- x$diss))
        out <- c(out, cl_validity.default(x, d))
    class(out) <- "cl_validity"
    out
}
## Package cluster: diana().
cl_validity.diana <-
function(x, ...)
{
    out <- list("Divisive coefficient" = x$dc)
    ## According to the docs, diana objects always have a diss
    ## component, but let's be defensive ...
    if(!is.null(d <- x$diss))
        out <- c(out, cl_validity.default(x, d))
    class(out) <- "cl_validity"
    out
}

## Package clue: (virtual) class "cl_partition".
cl_validity.cl_partition <-
function(x, ...)
    cl_validity(.get_representation(x), ...)
## Package clue: class pclust.
## So that this works for all classes extending pclust ...
cl_validity.pclust <-
function(x, ...)
    x$validity

print.cl_validity <-
function(x, ...)
{
    for(nm in names(x))
        cat(nm, ": ", x[[nm]], "\n", sep = "")
    invisible(x)
}

.cl_validity_partition_d_a_f <-
function(m, d)
{
    ## "Dissimilarity accounted for".
    ## Internal function for computing 1 - a / mean(d), where the
    ## "average within dissimilarity" a is given by
    ##   \frac{\sum_{i,j} \sum_k m_{ik}m_{jk} d(i,j)}
    ##        {\sum_{i,j} \sum_k m_{ik}m_{jk}}
    ## where m is the membership matrix and d a *symmetric* matrix of
    ## dissimilarities.
    within_sums <-
        rowSums(sapply(seq_len(ncol(m)),
                       function(k) {
                           z <- m[, k]
                           w <- outer(z, z, "*")
                           c(sum(w * d), sum(w))
                       }))
    average_within_d <- within_sums[1L] / within_sums[2L]
    1 - average_within_d / mean(d)
}

.cl_validity_hierarchy_variance_a_f <-
function(u, d)
{
    ## *Variance accounted for*.
    ## See e.g. Hubert, Arabie, & Meulman (2006), The structural
    ## representation of proximity matrices with MATLAB:
    ## variance_accounted_for = 
    ##   1 - \frac{\sum_{i < j} (d_{ij} - u_{ij}) ^ 2}
    ##            {\sum_{i < j} (d_{ij} - mean(d)) ^ 2}
    ## As this can be arbitrarily negative, we cut at 0.
    
    max(1 - sum((d - u) ^ 2) / sum((d - mean(d)) ^ 2), 0)
}

.cl_validity_hierarchy_deviance_a_f <-
function(u, d)
{
    ## *Deviance accounted for* (i.e., absolute deviation).
    ## See e.g. Smith (2001), Constructing ultrametric and additive
    ## trees based on the ${L}_1$ norm, Journal of Classification.
    ## deviance_accounted_for = 
    ##   1 - \frac{\sum_{i < j} |d_{ij} - u_{ij}|}
    ##            {\sum_{i < j} |d_{ij} - median(d)|}
    ## As this can be arbitrarily negative, we cut at 0.

    max(1 - sum(abs(d - u)) / sum(abs(d - median(d))), 0)
}

## Silhouette methods

silhouette.cl_partition <-
function(x, ...)
    silhouette(.get_representation(x), ...)

silhouette.cl_pclust <-
function(x, ...)
    x$silhouette

