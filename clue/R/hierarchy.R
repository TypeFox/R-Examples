### * is.cl_hierarchy

## Determine whether an object is a hierarchy.
## Note that hierarchies are n-trees, which can naturally be represented
## by their classes (as done via cl_classes()) or internal ultrametric
## obtained by assigning height one to all splits (as done by
## .cl_ultrametric_from_classes()). 
## We typically used the latter, but note that this is an *internal*
## reprsentation.
## User-level, cl_dendrogram objects are indexed hierarchies, and
## cl_hierarchy objects are n-trees.  The latter can be "converted" into
## the former (using height one splits) via as.cl_dendrogram().

is.cl_hierarchy <-
function(x)
    UseMethod("is.cl_hierarchy")

## Default method.
is.cl_hierarchy.default <- .false

## Package stats: hclust().
is.cl_hierarchy.hclust <-
function(x)
    !is.unsorted(x$height)

## Package cluster: agnes() and diana() give objects inheriting from
## class "twins".
is.cl_hierarchy.twins <- .true
## Package cluster: mona().
is.cl_hierarchy.mona <- .true

## Package ape: class "phylo".
is.cl_hierarchy.phylo <-
function(x)
    ape::is.ultrametric(x)

## Package clue: (virtual) class "cl_hierarchy".
## Note that "raw" cl_ultrametric objects are *not* hierarchies, as
## these are meant for numeric computations.
## <FIXME>
## Is this really a good idea?
## We can as.hclust() a cl_dendrogram and then it is a cl_hierarchy ...
## </FIXME>
is.cl_hierarchy.cl_hierarchy <- .true

### * as.cl_hierarchy

## Note that cl_hierarchy conceptually is a virtual class, so there are
## no prototypes and no cl_hierarchy() creator.

.cl_hierarchy_classes <- "cl_hierarchy"

as.cl_hierarchy <-
function(x)
{
    if(is.cl_hierarchy(x)) {
        if(!inherits(x, "cl_hierarchy"))
            .make_container(x, .cl_hierarchy_classes)
        else
            x
    }
    else
        .make_container(as.cl_ultrametric(x),
                        .cl_hierarchy_classes)
}

### * print.cl_hierarchy

print.cl_hierarchy <-
function(x, ...)
    .print_container(x, "cl_hierarchy", ...)

### * Complex.cl_hierarchy

## No Complex() for any kind of hierarchy.
Complex.cl_hierarchy <-
function(z)
    stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                  .Generic, .Class),
         domain = NA)

### * Math.cl_hierarchy

## No Math() for any kind of hierarchy.
Math.cl_hierarchy <-
function(x, ...)
    stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                  .Generic, .Class),
         domain = NA)

### * Ops.cl_hierarchy

Ops.cl_hierarchy <-
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

    if(n_of_objects(e1) != n_of_objects(e2))
        stop("Hierarchies must have the same number of objects.")

    c1 <- cl_classes(e1)
    c2 <- cl_classes(e2)

    switch(.Generic,
           "<=" = all(is.finite(match(c1, c2))),
           "<"  = all(is.finite(match(c1, c2))) && any(is.na(match(c2, c1))),
           ">=" = all(is.finite(match(c2, c1))),
           ">"  = all(is.finite(match(c2, c1))) && any(is.na(match(c1, c2))),
           "==" = all(is.finite(match(c1, c2))) && all(is.finite(match(c2, c1))),
           "!=" = any(is.na(match(c1, c2))) || any(is.na(match(c2, c1))))
}

### * Summary.cl_hierarchy

## <NOTE>
## This is really the same as Summary.cl_partition().
## </NOTE>

Summary.cl_hierarchy <-
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

### * as.hclust.cl_hierarchy

as.hclust.cl_hierarchy <-
function(x, ...)
    as.hclust(.get_representation(x), ...)

### * is.cl_dendrogram

## <FIXME>
## Once we have cl_dendrogram testing, we can simplify cl_hierarchy
## testing.  E.g.,
##   is.cl_hierachy.default <- is.cl_dendrogram
## should be ok, and we can add cl_hierarchy predicates for hierarchies
## which are not dendrograms on top of that.
## </FIXME>

is.cl_dendrogram <-
function(x)
    UseMethod("is.cl_dendrogram")
## Default method.
is.cl_dendrogram.default <- .false
## Package stats: hclust().
is.cl_dendrogram.hclust <-
function(x)
    !is.unsorted(x$height)
## Package cluster: agnes() and diana() give objects inheriting from
## class "twins".
is.cl_dendrogram.twins <- .true
## Package cluster: mona().
is.cl_dendrogram.mona <- .true
## Package ape: class "phylo".
is.cl_dendrogram.phylo <-
function(x)
    ape::is.ultrametric(x)
## (We could also support ape's class "matching" via coercion to class
## "phylo".)
## Package clue: (virtual) class "cl_dendrogram".
is.cl_dendrogram.cl_dendrogram <- .true

### * as.cl_dendrogram

.cl_dendrogram_classes <- c("cl_dendrogram", "cl_hierarchy")

as.cl_dendrogram <-
function(x)
{
    if(is.cl_dendrogram(x)) {
        if(!inherits(x, "cl_dendrogram"))
            .make_container(x, .cl_dendrogram_classes)
        else
            x
    }
    else
        .make_container(as.cl_ultrametric(x),
                        .cl_dendrogram_classes)
}

### * print.cl_dendrogram

print.cl_dendrogram <-
function(x, ...)
    .print_container(x, "cl_dendrogram", ...)

### * plot.cl_dendrogram

plot.cl_dendrogram <-
function(x, ...)
    plot(cl_ultrametric(.get_representation(x)), ...)

### * Group methods for cl_dendrogram objects.

Ops.cl_dendrogram <-
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

    u1 <- cl_ultrametric(e1)
    u2 <- cl_ultrametric(e2)
    if(length(u1) != length(u2))
        stop("Dendrograms must have the same number of objects.")
    switch(.Generic,
           "<=" = all(u1 <= u2),
           "<"  = all(u1 <= u2) && any(u1 < u2),
           ">=" = all(u1 >= u2),
           ">"  = all(u1 >= u2) && any(u1 > u2),
           "==" = all(u1 == u2),
           "!=" = any(u1 != u2))
}

### * Summary.cl_dendrogram

## <NOTE>
## This is really the same as Summary.cl_hierarchy() ...
## We cannot really call the poset specific internal meet and join
## functions from here as e.g. max(D, H) (D a dendrogram, H an n-tree)
## should use the n-tree poset functions ...
## However, dispatch for cl_dendrogram should not be needed if we also
## dispatch on cl_hierarchy ...
## </NOTE>

## Summary.cl_dendrogram <-
## function(..., na.rm = FALSE)
## {
##     ok <- switch(.Generic, max = , min = , range = TRUE, FALSE)
##     if(!ok)
##         stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
##                       .Generic, .Class))
##     args <- list(...)
##     switch(.Generic,
##            "min" = cl_meet(cl_ensemble(list = args)),
##            "max" = cl_join(cl_ensemble(list = args)),
##            "range" = {
##                cl_ensemble(min = cl_meet(cl_ensemble(list = args)),
##                            max = cl_join(cl_ensemble(list = args)))
##            })
## }

### * as.hclust.cl_dendrogram

## <NOTE>
## This is really the same as as.hclust.cl_hierarchy() ...
## Dispatch for cl_dendrogram should not be needed if we also dispatch
## on cl_hierarchy ...
## </NOTE>

## as.hclust.cl_dendrogram <-
## function(x, ...)
##     as.hclust(.get_representation(x), ...)

### ** cut.cl_dendrogram

## Not perfect as this perhaps return something more "classed" in the
## spirit of clue ...

cut.cl_dendrogram <-
function(x, ...)
    cutree(as.hclust(x), ...)

### * Utilities

## To turn a mona object into a cl_dendrogram, we need to be able to
## compute its associated ultrametric.  Hence, provide a cophenetic()
## method for mona objects ...
cophenetic.mona <-
function(x)
{
    no <- length(x$order)
    ns <- max(x$step) + 1

    m <- matrix(NA, no, no)
    FOO <- function(ind, step, s) {
        if(length(ind) <= 1) return()
        grp <- c(0, cumsum(step == s))
        ind <- split(ind, grp)
        len <- length(ind)
        for(a in seq_len(len)) {
            for(b in seq(from = 1, length.out = a - 1)) {
                ## Need both as we currently cannot assume that the
                ## indices are sorted.  Alternatively, work with the
                ## sequence from one to the number of objects, and
                ## reorder at the end ...
                m[ind[[a]], ind[[b]]] <<- s
                m[ind[[b]], ind[[a]]] <<- s
            }
        }
        ind <- ind[sapply(ind, length) > 1]
        pos <- which(step == s)
        step <- split(step[-pos], grp[-1][-pos])
        if(is.null(step)) return()
        for(a in seq_along(ind))
            FOO(ind[[a]], step[[a]], s + 1)
    }
    
    FOO(x$order, x$step, 1)
    m[is.na(m)] <- ns
    m <- ns - m
    rownames(m) <- rownames(x$data)
    as.dist(m)
}

## And while we're at it ...
## (Of course, as.hclust() should really "know" that a cophenetic()
## method is available ...)
as.hclust.mona <-
function(x, ...)
    hclust(cophenetic(x), "single")
    

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
