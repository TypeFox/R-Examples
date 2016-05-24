## <FIXME>
## Maybe add support for "auto" type (class_ids when predicting from a
## hard, memberships when predicting from a soft partition) eventually.
## </FIXME>

cl_predict <-
function(object, newdata = NULL,
         type = c("class_ids", "memberships"), ...)
    UseMethod("cl_predict")

## Default method.
## Should also work for kcca() from package flexclust.
cl_predict.default <-
function(object, newdata = NULL,
         type = c("class_ids", "memberships"), ...)
    .as_cl_class_ids_or_membership(predict(object, newdata, ...), type)

## Package stats: kmeans() (R 2.1.0 or better).
cl_predict.kmeans <-
function(object, newdata = NULL,
         type = c("class_ids", "memberships"), ...)
{
    if(is.null(newdata))
        return(.cl_class_ids_or_membership(object, type))
    d <- .rxdist(newdata, object$centers)
    .as_cl_class_ids_or_membership(max.col(-d), type)
}

## Package cluster:
## * fanny() cannot make "new" predictions.
## * clara() gives medoids, and takes metric data using Euclidean or
##   Manhattan dissimilarities (and we can figure out which by looking
##   at the call and the default values).
## * pam() gives medoids, but might have been called with dissimilarity
##   data, so is tricky.  We can always find out which by looking at the
##   medoids: as in the dissimilarity input case this is a vector of
##   class labels, and a matrix with in each row the coordinates of one
##   medoid otherwise.  We then still need to figure out whether
##   Euclidean or Manhattan distances were used by looking at the call
##   and the default values.
## Both pam() and clara() show that the interfaces could be improved to
## accomodate modern needs, e.g., for bagging.

cl_predict.fanny <-
function(object, newdata = NULL,
         type = c("class_ids", "memberships"), ...)
{
    if(is.null(newdata))
        return(.cl_class_ids_or_membership(object, type))
    stop("Cannot make new predictions.")
}

cl_predict.clara <-
function(object, newdata = NULL,
         type = c("class_ids", "memberships"), ...)
{
    if(is.null(newdata))
        return(.cl_class_ids_or_membership(object, type))
    ## <FIXME>
    ## Add support eventually ...
    if(identical(object$call$stand, TRUE))
        warning("Standardization is currently not supported.")
    ## </FIXME>
    method <- object$call$metric
    if(is.null(method)) {
        ## Not given in the call, hence use default value.
        method <- formals(cluster::clara)$metric
        ## (Or hard-wire the default value: "euclidean".)
    }
    d <- .rxdist(newdata, object$medoids, method)
    .as_cl_class_ids_or_membership(max.col(-d), type)
}

cl_predict.pam <-
function(object, newdata = NULL,
         type = c("class_ids", "memberships"), ...)
{
    if(is.null(newdata))
        return(.cl_class_ids_or_membership(object, type))
    prototypes <- object$medoids
    if(!is.matrix(prototypes))
        stop("Cannot make new predictions.")
    ## <FIXME>
    ## Add support eventually ...
    if(identical(object$call$stand, TRUE))
        warning("Standardization is currently not supported.")
    ## </FIXME>
    method <- object$call$metric
    if(is.null(method)) {
        ## Not given in the call, hence use default value.
        method <- formals(cluster::pam)$metric
        ## (Or hard-wire the default value: "euclidean".)
    }
    d <- .rxdist(newdata, object$medoids, method)
    .as_cl_class_ids_or_membership(max.col(-d), type)
}

## Package RWeka: clusterers return objects inheriting from
## "Weka_clusterer".
cl_predict.Weka_clusterer <-
function(object, newdata = NULL,
         type = c("class_ids", "memberships"), ...)
{
    if(is.null(newdata))
        return(.cl_class_ids_or_membership(object, type))
    .as_cl_class_ids_or_membership(predict(object, newdata = newdata,
                                           type = type, ...),
                                   type)
}

## Package cba: ccfkms().
cl_predict.ccfkms <-
function(object, newdata = NULL,
         type = c("class_ids", "memberships"), ...)
{
    if(is.null(newdata))
        return(.cl_class_ids_or_membership(object, type))
    .as_cl_class_ids_or_membership(as.vector(predict(object,
                                                     newdata)$cl),
                                   type)
}
## Package cba: rockCluster() returns objects of class "rock".
## If x is a Rock object, fitted(x) and predict(x, newdata) can result
## in missing classifications, as
##   In the case a 'drop' value greater than zero is specified, all
##   clusters with size equal or less than this value are removed from
##   the classifier.  Especially, 'fitted' uses a threshold of one
##   because for singleton clusters the neighborhood is empty.
cl_predict.rock <-
function(object, newdata = NULL,
         type = c("class_ids", "memberships"), ...)
{
    if(is.null(newdata)) newdata <- object$x
    ids <- as.vector(predict(object, newdata, ...)$cl)
    .as_cl_class_ids_or_membership(ids, type)
}

## Package cclust: cclust().
cl_predict.cclust <-
function(object, newdata = NULL,
         type = c("class_ids", "memberships"), ...)
{
    ## Package cclust provides predict.cclust() which returns (again) an
    ## object of class "cclust", but does not give the labels of the
    ## original data in case no new data are given.
    if(is.null(newdata))
        return(.cl_class_ids_or_membership(object, type))
    .as_cl_class_ids_or_membership(predict(object, newdata), type)
}

## Package e1071: cmeans() gives objects of class "fclust".
cl_predict.fclust <-
function(object, newdata = NULL,
         type = c("class_ids", "memberships"), ...)
{
    if(is.null(newdata))
        return(.cl_class_ids_or_membership(object, type))

    ## Note that the 'fclust' objects returned by cmeans() do not always
    ## directly contain the information on the fuzzification parameter m
    ## and the distance (Euclidean/Manhattan) employed, so we have to
    ## engineer this from the matched call and the default arguments.
    nms <- names(object$call)
    ## Note that we cannot directly use object$call$m, as this could
    ## give the 'method' argument if 'm' was not given.
    m <- if("m" %in% nms)
        object$call$m
    else {
        ## Not given in the call, hence use default value.
        formals(e1071::cmeans)$m
        ## (Or hard-wire the default value: 2.)
    }
    method <- if("dist" %in% nms)
        object$call$dist
    else {
        ## Not given in the call, hence use default value.
        formals(e1071::cmeans)$dist
        ## (Or hard-wire the default value: "euclidean".)
    }

    d <- .rxdist(newdata, object$centers, method)
    power <- c(m, if(method == "euclidean") 2 else 1)
    M <- .memberships_from_cross_dissimilarities(d, power)
    .as_cl_class_ids_or_membership(M, type)
}

## Package e1071: cshell().
cl_predict.cshell <-
function(object, newdata = NULL,
         type = c("class_ids", "memberships"), ...)
{
    if(is.null(newdata))
        return(.cl_class_ids_or_membership(object, type))

    ## Not surprisingly, this is rather similar to what we do for fclust
    ## objects.  Only dissimiliraties (and exponents) need to be
    ## computed differently ...
    nms <- names(object$call)
    m <- if("m" %in% nms)
        object$call$m
    else {
        ## Not given in the call, hence use default value.
        formals(e1071::cshell)$m
        ## (Or hard-wire the default value: 2.)
    }
    method <- if("dist" %in% nms)
        object$call$dist
    else {
        ## Not given in the call, hence use default value.
        formals(e1071::cshell)$dist
        ## (Or hard-wire the default value: "euclidean".)
    }

    d <- .rxdist(newdata, object$centers, method)
    d <- sweep(d, 2, object$radius) ^ 2
    M <- .memberships_from_cross_dissimilarities(d, m)
    .as_cl_class_ids_or_membership(M, type)
}

## Package e1071: bclust().
## <NOTE>
## One might argue that it would be better to use the 'dist.method'
## employed for the hierarchical clustering, but it seems that class
## labels ("clusters") are always assigned using Euclidean distances.
cl_predict.bclust <- cl_predict.kmeans
## </NOTE>

## Package flexclust: kcca() returns objects of S4 class "kcca" which
## extends S4 class "flexclust".
cl_predict.kcca <- cl_predict.default

## Package flexmix: class "flexmix".
cl_predict.flexmix <-
function(object, newdata = NULL,
         type = c("class_ids", "memberships"), ...)
{
    if(is.null(newdata))
        return(.cl_class_ids_or_membership(object, type))
    .as_cl_class_ids_or_membership(modeltools::posterior(object,
                                                         newdata,
                                                         ...),
                                   type)
}

## Package mclust: Mclust().
cl_predict.Mclust <-
function(object, newdata = NULL,
         type = c("class_ids", "memberships"), ...)
{
    if(is.null(newdata))
        return(.cl_class_ids_or_membership(object, type))

    pred <- predict(object, newdata, ...)
    type <- match.arg(type)
    if(type == "class_ids")
        as.cl_class_ids(pred$classification)
    else
        as.cl_membership(pred$z)
}


## Package movMF: movMF().
cl_predict.movMF <- cl_predict.Weka_clusterer

## Package clue: pclust().
cl_predict.pclust <-
function(object, newdata = NULL,
         type = c("class_ids", "memberships"), ...)
{
    if(is.null(newdata))
        return(.cl_class_ids_or_membership(object, type))

    d <- object$family$D(newdata, object$prototypes)
    power <- c(object$m, object$family$e)
    M <- .memberships_from_cross_dissimilarities(d, power)
    .as_cl_class_ids_or_membership(M, type)
}

## Package clue: (virtual) class "cl_partition".
cl_predict.cl_partition <-
function(object, newdata = NULL,
         type = c("class_ids", "memberships"), ...)
    cl_predict(.get_representation(object), newdata = newdata, type, ...)

## Internal helpers: this looks a bit silly, but makes the rest of the
## code look nicer ...

.cl_class_ids_or_membership <-
function(x, type = c("class_ids", "memberships"))
{
    type <- match.arg(type)

    if(type == "class_ids")
        cl_class_ids(x)
    else
        cl_membership(x)
}

.as_cl_class_ids_or_membership <-
function(x, type = c("class_ids", "memberships"))
{
    type <- match.arg(type)

    if(type == "class_ids") {
        if(is.matrix(x)) {
            ## Same as for cl_class_ids.cl_membership().
            as.cl_class_ids(.structure(max.col(x), names = rownames(x)))
        }
        else
            as.cl_class_ids(x)
    }
    else
        as.cl_membership(x)
}
