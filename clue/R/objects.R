### * n_of_objects

## Get the number of objects in a clustering.

n_of_objects <-
function(x)
    UseMethod("n_of_objects")

### ** Default method.

n_of_objects.default <-
function(x)
    length(cl_class_ids(x))
## (Note that prior to R 2.1.0, kmeans() returned unclassed results,
## hence the best we can do for the *default* method is to look at a
## possibly existing "cluster" component.  Using the class ids incurs
## another round of method dispatch, but avoids code duplication.)

### ** Partitioning methods.

## Package stats: kmeans() (R 2.1.0 or better).
n_of_objects.kmeans <- n_of_objects.default

## Package cluster: clara(), fanny(), and pam() give objects of the
## respective class inheriting from class "partition".
n_of_objects.partition <- n_of_objects.default

## Package cclust: cclust().
n_of_objects.cclust <- n_of_objects.default

## Package e1071: cmeans() gives objects of class "fclust".
n_of_objects.fclust <-
function(x)
    nrow(x$membership)
## Package e1071: cshell().
n_of_objects.cshell <- n_of_objects.fclust
## Package e1071: bclust().
n_of_objects.bclust <- n_of_objects.default

## Package mclust: Mclust().
n_of_objects.Mclust <- n_of_objects.default

### ** Hierarchical methods.

## Package stats: hclust().
n_of_objects.hclust <-
function(x)
    length(x$order)

## Package cluster: agnes() and diana() give objects inheriting from
## class "twins".
n_of_objects.twins <- n_of_objects.hclust
## Package cluster: mona().
n_of_objects.mona <- n_of_objects.hclust

## Package ape: class "phylo".
n_of_objects.phylo <-
function(x)
    length(x$tip.label)

### ** Others.

## Package stats: class "dist".
n_of_objects.dist <-
function(x)
    attr(x, "Size")

## Package clue: Ensembles.
n_of_objects.cl_ensemble <-
function(x)
    attr(x, "n_of_objects")
## Package clue: Memberships.
n_of_objects.cl_membership <- nrow
## Package clue: pclust().
n_of_objects.pclust <- n_of_objects.default
## Package clue: Ultrametrics.
n_of_objects.cl_ultrametric <- n_of_objects.dist
## Package clue: (virtual) class "cl_partition".
n_of_objects.cl_partition <-
function(x)
    .get_property_from_object_or_representation(x, "n_of_objects")
## Package clue: (virtual) class "cl_hierarchy".
n_of_objects.cl_hierarchy <-
function(x)
    .get_property_from_object_or_representation(x, "n_of_objects")

### * cl_object_names

## Determine the names of the objects in a clustering if available; give
## NULL otherwise.  This is in sync with e.g. names() or dimnames(); au
## contraire, cl_object_labels() always gives labels even if no names
## are available.

cl_object_names <- 
function(x)
    UseMethod("cl_object_names")

## ** Default method.

cl_object_names.default <- function(x) names(cl_class_ids(x))

## ** Partitions.

## There is really nothing special we can currently do.
## Most partitioning functions return no information on object names.
## This includes classes
##   stats:      kmeans
##   cba:        ccfkms, rock
##   cclust:     cclust
##   e1071:      bclust
##   flexclust:  kcca
##   kernlab:    specc
##   mclust:     Mclust
## The algorithms for which things "work" all give named class ids.
##   RWeka:      Weka_clusterer
##   cluster:    clara fanny pam
##   e1071:      cclust cshell

## ** Hierarchies.

## Package stats: hclust().
cl_object_names.hclust <- function(x) x$labels

## Package cluster: agnes(), diana() and mona() all return an object
## which has an 'order.lab' component iff "the original observations
## were labelled".  We can use this together the the 'order' component
## to recreate the labels in their original order.  Note that we cannot
## rely on dissimilarity or data components being available.
cl_object_names.twins <-
function(x)
{
    if(!is.null(x$order.lab)) {
        out <- character(length = n_of_objects(x))
        out[x$order] <- x$order.lab
        out
    }
    else
        NULL
}
cl_object_names.mona <- cl_object_names.twins

## Package ape: class "phylo".
cl_object_names.phylo <-
function(x)
    x$tip.label

## ** Others.

## Package stats: class "dist".
## (Raw object dissimilarities.)
cl_object_names.dist <-
function(x)
    attr(x, "Labels")

## Package clue: memberships.
cl_object_names.cl_membership <-
function(x)
    rownames(x)
## Package clue: ultrametrics.
cl_object_names.cl_ultrametric <-
function(x)
    attr(x, "Labels")
## Package clue: (virtual) class "cl_partition".
cl_object_names.cl_partition <-
function(x)
    cl_object_names(.get_representation(x))
## Package clue: (virtual) class "cl_hierarchy".
cl_object_names.cl_hierarchy <-
function(x)
    cl_object_names(.get_representation(x))
## Package clue: ensembles.
cl_object_names.cl_ensemble <-
function(x)
{
    nms <- lapply(x, cl_object_names)
    ind <- which(sapply(nms, length) > 0)
    if(any(ind)) nms[[ind[1L]]] else NULL
}

### * cl_object_labels

cl_object_labels <-
function(x)
{
    if(is.null(out <- cl_object_names(x)))
        out <- as.character(seq_len(n_of_objects(x)))
    out
}

### * cl_object_dissimilarities

## Extract object dissimilarities from R objects containing such: this
## includes objects directly inheriting from "dist" as well as
## dendrograms or additive trees.

cl_object_dissimilarities <-
function(x)
{
    ## Keep this in sync with .has_object_dissimilarities().
    if(is.cl_dendrogram(x))
        cl_ultrametric(x)
    else if(inherits(x, "dist"))
        x
    else
        stop("Cannot extract object dissimilarities")
}

.has_object_dissimilarities <-
function(x)
{
    ## Keep this in sync with cl_object_dissimilarities().
    is.cl_dendrogram(x) || inherits(x, "dist")
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***

