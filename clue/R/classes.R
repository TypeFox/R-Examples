cl_classes <-
function(x)
    UseMethod("cl_classes")

cl_classes.default <-
function(x)
{
    ## Be nice to users ...
    if(is.cl_partition(x))
        cl_classes(as.cl_partition(x))
    else if(is.cl_dendrogram(x))
        cl_classes(as.cl_dendrogram(x))
    else
        stop("Can only determine classes of partitions or hierarchies.")
}

cl_classes.cl_partition <-
function(x)
{
    n <- n_of_objects(x)
    out <- split(seq_len(n), cl_class_ids(x))
    class(out) <- c("cl_classes_of_partition_of_objects",
                    "cl_classes_of_objects")
    attr(out, "n_of_objects") <- n
    attr(out, "labels") <- cl_object_labels(x)
    out
}

cl_classes.cl_hierarchy <-
function(x)
{
    ## Assume a valid hierarchy/dendrogram.
    x <- as.hclust(x)
    n <- n_of_objects(x)
    labels <- seq_len(n)
    ## Only use the "maximal" partitions for each height (relevant in
    ## case of non-binary trees).
    groups <- cutree(x, h = unique(c(0, x$height)))
    ## Give a list with the (unique) sets of numbers of the objects.
    out <- unique(unlist(sapply(split(groups, col(groups)),
                                function(k) split(labels, k)),
                         recursive = FALSE,
                         use.names = FALSE))
    ## Preserve labels if possible, and re-order according to
    ## cardinality.
    out <- out[order(sapply(out, length))]
    class(out) <- c("cl_classes_of_hierarchy_of_objects",
                    "cl_classes_of_objects")
    attr(out, "n_of_objects") <- n
    attr(out, "labels") <- cl_object_labels(x)
    out
}

## Be nice to users of ultrametric fitters ... which should really fit
## dendrograms (which inherit from hierarchies).
cl_classes.cl_ultrametric <- cl_classes.cl_hierarchy

print.cl_classes_of_partition_of_objects <-
function(x, ...)
{
    labels <- attr(x, "labels")
    y <- lapply(x, function(i) paste(labels[i], collapse = ", "))
    writeLines(formatDL(names(x), sprintf("{%s}", unlist(y)),
                        style = "list", ...))
    invisible(x)
}

print.cl_classes_of_hierarchy_of_objects <-
function(x, ...)
{
    labels <- attr(x, "labels")
    y <- lapply(x, function(i) paste(labels[i], collapse = ", "))
    y <- strwrap(sprintf("{%s},", unlist(y)), exdent = 2)
    y[length(y)] <- sub(",$", "", y[length(y)])
    writeLines(y)
    invisible(x)
}
