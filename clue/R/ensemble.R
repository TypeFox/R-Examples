cl_ensemble <-
function(..., list = NULL)
{
    clusterings <- c(list(...), list)

    if(!length(clusterings)) {
        ## Return an empty cl_ensemble.
        ## In this case, we cannot additionally know whether it contains
        ## partitions or hierarchies ...
        attr(clusterings, "n_of_objects") <- as.integer(NA)
        class(clusterings) <- "cl_ensemble"
        return(clusterings)
    }

    ## Previously, we used to require that the elements of the ensemble
    ## either all be partitions, or all be hierarchies.  We no longer do
    ## this, as it makes sense to also allow e.g. object dissimilarities
    ## (raw "dist" objects or additive distances) as elements (e.g.,
    ## when computing proximities), and it is rather cumbersome to
    ## decide in advance which combinations of elements might be useful
    ## and hence should be allowed.  All we enforce is that all elements
    ## correspond to the same number of objects (as we typically cannot
    ## verify that they relate to the *same* objects).  For "pure"
    ## ensembles of partitions or hierarchies we add additional class
    ## information.

    if(all(sapply(clusterings, is.cl_partition)))
        class(clusterings) <- c("cl_partition_ensemble", "cl_ensemble")
    else if(all(sapply(clusterings, is.cl_dendrogram)))
        class(clusterings) <- c("cl_dendrogram_ensemble",
                                "cl_hierarchy_ensemble", "cl_ensemble")
    else if(all(sapply(clusterings, is.cl_hierarchy)))
        class(clusterings) <- c("cl_hierarchy_ensemble", "cl_ensemble")
    else
        class(clusterings) <- "cl_ensemble"

    n <- sapply(clusterings, n_of_objects)
    if(any(diff(n)))
        stop("All elements must have the same number of objects.")
    attr(clusterings, "n_of_objects") <- as.integer(n[1L])

    clusterings
}

is.cl_ensemble <-
function(x)
    inherits(x, "cl_ensemble")

## <NOTE>
## In the old days, kmeans() results were unclassed lists, hence such
## objects were taken as representing a single clustering.  Nowadays, we
## take these as lists of clusterings.
as.cl_ensemble <-
function(x)
{
    if(is.cl_ensemble(x)) x
    else if(is.list(x) && !is.object(x)) cl_ensemble(list = x)
    else cl_ensemble(x)
}
## </NOTE>

c.cl_ensemble <-
function(..., recursive = FALSE)
{
    clusterings <- unlist(lapply(list(...), as.cl_ensemble),
                          recursive = FALSE)
    cl_ensemble(list = clusterings)
}

"[.cl_ensemble" <-
function(x, i)
{
    ## Make subscripting empty ensembles a noop.
    if(length(x) == 0L) return(x)
    cl_ensemble(list = NextMethod("["))
}

rep.cl_ensemble <-
function(x, times, ...)
    cl_ensemble(list = NextMethod("rep"))

print.cl_partition_ensemble <-
function(x, ...)
{
    msg <- sprintf(ngettext(length(x),
                            "An ensemble of %d partition of %d objects.",
                            "An ensemble of %d partitions of %d objects."),
                   length(x), n_of_objects(x))
    writeLines(strwrap(msg))
    invisible(x)
}

Summary.cl_partition_ensemble <-
function(..., na.rm = FALSE)
{
    ok <- switch(.Generic, max = , min = , range = TRUE, FALSE)
    if(!ok)
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class),
             domain = NA)
    args <- list(...)
    ## Combine the given partition ensembles.
    x <- do.call(c, args)
    switch(.Generic,
           "min" = cl_meet(x),
           "max" = cl_join(x),
           "range" = cl_ensemble(min = cl_meet(x), max = cl_join(x)))
}

print.cl_dendrogram_ensemble <-
function(x, ...)
{
    msg <- sprintf(ngettext(length(x),
                            "An ensemble of %d dendrogram of %d objects.",
                            "An ensemble of %d dendrograms of %d objects."),
                   length(x), n_of_objects(x))
    writeLines(strwrap(msg))
    invisible(x)
}

print.cl_hierarchy_ensemble <-
function(x, ...)
{
    msg <- sprintf(ngettext(length(x),
                            "An ensemble of %d hierarchy of %d objects.",
                            "An ensemble of %d hierarchies of %d objects."),
                   length(x), n_of_objects(x))
    writeLines(strwrap(msg))
    invisible(x)
}

print.cl_ensemble <-
function(x, ...)
{
    writeLines(sprintf(ngettext(length(x),
                                "An ensemble with %d element.",
                                "An ensemble with %d elements."),
                       length(x)))
    invisible(x)
}

plot.cl_ensemble <-
function(x, ..., main = NULL, layout = NULL)
{
    if(!is.cl_ensemble(x))
        stop("Wrong class.")

    ## What we can definitely plot is are cl_addtree, cl_dendrogram and
    ## cl_ultrametric objects.  (We could also add simple methods for
    ## plotting raw dissimilarities, but of course seriation::dissplot()
    ## would be the thing to use.)  What we cannot reasonably plot is
    ## partitions (in particular, as these do not know about the
    ## underlying dissimilarities.  But then we could perhaps provide
    ## silhoutte plots etc for ensembles of partitions ...
    ## <FIXME>
    ## Think about this.
    ## </FIXME>
    ## So let us check for the things we can plot.

    ## (Note that currently there is neither is.cl_ultrametric() nor
    ## is.cl_addtree().)
    ok <- sapply(x,
                 function(e)
                 (is.cl_dendrogram(e) ||
                  inherits(e, c("cl_addtree", "cl_ultrametric"))))
    if(!all(ok))
        stop(gettextf("Plotting not available for elements %s of the ensemble.",
                      paste(which(!ok), collapse = " ")),
             domain = NA)

    ## Prefer dendrogram plot methods to those for hclust objects.
    ind <- sapply(x, is.cl_dendrogram)
    if(any(ind))
        x[ind] <- lapply(x, as.cl_dendrogram)

    ## Now the usual layouting ... same as for plotting relation
    ## ensembles.

    ## Number of elements.
    n <- length(x)
    ## Layout.
    byrow <- TRUE
    if(is.null(layout)) {
        nc <- ceiling(sqrt(n))
        nr <- ceiling(n / nc)
    }
    else {
        layout <- c(as.list(layout), byrow)[seq_len(3)]
        if(is.null(names(layout)))
            names(layout) <- c("nr", "nc", "byrow")
        nr <- layout[["nr"]]
        nc <- layout[["nc"]]
        byrow <- layout[["byrow"]]
    }
    op <- if(byrow)
        par(mfrow = c(nr, nc))
    else
        par(mfcol = c(nr, nc))
    on.exit(par(op))

    ## Try recycling main (might want the same for others as well).
    if(!is.list(main)) {
        main <- if(is.null(main))
            vector("list", length = n)
        else
            rep.int(as.list(main), n)
    }
    
    for(i in seq_along(x))
        plot(x[[i]], main = main[[i]], ...)
}
                   
unique.cl_ensemble <-
function(x, incomparables = FALSE, ...)
    cl_ensemble(list = NextMethod("unique"))

.cl_ensemble_type <-
function(x)
{
    if(inherits(x, "cl_partition_ensemble"))
        "partition"
    else if(inherits(x, "cl_hierarchy_ensemble"))
        "hierarchy"
    else
        NULL
}
