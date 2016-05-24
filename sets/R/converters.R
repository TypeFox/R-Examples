### * converters

### SET CONVERTERS
### We basically have make_set_with_order converters
### which return the resulting set and the original order
### to allow the permutation of some associated meta information
### according to the new order.
### (Examples: memberships of generalized sets and incidences of relations.)
### The as.set converter calls these internally and
### strips the ordering information.

make_set_with_order <-
function(x)
    UseMethod("make_set_with_order")

make_set_with_order.default <-
function(x)
    stop("Not implemented.")

make_set_with_order.NULL <-
function(x)
    make_set_with_order(list())

make_set_with_order.set <-
function(x)
    .make_set_with_order(x)

make_set_with_order.gset <-
function(x)
{
    attr(x, "memberships") <- NULL
    .make_set_with_order(.make_set_from_list(as.list(.get_support(x))))
}

make_set_with_order.cset <-
function(x)
    make_set_with_order(as.gset(x))

make_set_with_order.numeric <-
make_set_with_order.integer <-
function(x)
{
    x <- unique(x)
    O <- order(x)
    .make_set_with_order(.make_set_from_list(as.list(x[O])), O)
}

make_set_with_order.factor <-
make_set_with_order.character <-
make_set_with_order.logical <-
make_set_with_order.ordered <-
make_set_with_order.tuple <-
make_set_with_order.list <-
function(x)
{
    x <- .list_unique(x)
    O <- .list_order(x)
    .make_set_with_order(.make_set_from_list(x[O]), O)
}

make_set_with_order.matrix <-
make_set_with_order.data.frame <-
function(x) {
    x <- unique(x)
    O <- do.call(order, x)
    n <- rownames(x)
    if (is.null(n)) n <- seq_len(nrow(x))
    .make_set_with_order(.make_set_from_list(lapply(split(x, n), as.tuple)), O)
}

.make_set_with_order <-
function(set, order = seq_along(set))
    list(set = set, order = order)

### High-Level converter

as.set <-
function(x)
    make_set_with_order(x)$set

### canonicalizer

canonicalize_set_and_mapping <-
function(x, mapping = NULL, margin = NULL)
{
    x <- make_set_with_order(x)
    if (!is.null(mapping))
        mapping <- if (is.array(mapping) || is.data.frame(mapping)) {
            D <- dim(mapping)
            L <- length(x$set)
            if (is.null(margin))
                margin <- which(D == L)
            permute <- rep.int(list(seq_len(L)), length(D))
            permute[margin] <- rep.int(list(x$order), length(margin))
            do.call("[", c(list(mapping), permute, list(drop = FALSE)))
        } else
            mapping[x$order]

    list(set = x$set, mapping = mapping, order = x$order)
}

###

as.list.set <-
function(x, ...)
    .as.list(x)

### gset converters

as.gset <-
function(x)
    UseMethod("as.gset")

as.gset.default <-
function(x)
    gset(x)

as.gset.cset <-
function(x)
{
    attr(x, "orderfun") <- NULL
    attr(x, "matchfun") <- NULL
    class(x) <- if (!is.null(attr(x, "memberships")))
        c("gset")
    else
        c("set", "gset")
    x
}

as.gset.gset <- identity

as.gset.tuple <-
function(x)
    as.gset(as.list(x))

as.gset.numeric <- function(x) .as.gset.atomic(x, as.numeric)
as.gset.character <- function(x) .as.gset.atomic(x, as.character)
as.gset.factor <- function(x) .as.gset.atomic(x, as.factor)
as.gset.ordered <- function(x) .as.gset.atomic(x, as.ordered)
as.gset.integer <- function(x) .as.gset.atomic(x, as.integer)
as.gset.logical <- function(x) .as.gset.atomic(x, as.logical)

.as.gset.atomic <-
function(x, FUN)
{
    tab <- table(x)
    .make_gset_from_support_and_memberships(FUN(c(names(tab), NA)),
                                            c(as.vector(tab), sum(is.na(x))))
}


as.gset.list <-
function(x)
{
    uni <- .list_sort(.list_unique(x))
    count <- table(.exact_match(x, uni))
    .make_gset_from_support_and_memberships(uni, count)
}

as.gset.data.frame <-
as.gset.matrix <-
function(x)
{
    n <- rownames(x)
    if (is.null(n)) n <- seq_len(nrow(x))
    as.gset(lapply(split(x, n), as.tuple))
}

as.list.gset <-
function(x, ...)
    .as.list(x)

### tuple converters

as.tuple <-
function(x)
    UseMethod("as.tuple")

as.tuple.default <-
function(x)
    stop("Not implemented!")

as.tuple.tuple <- identity

as.tuple.numeric <-
as.tuple.factor <-
as.tuple.character <-
as.tuple.integer <-
as.tuple.ordered <-
as.tuple.logical <-
function(x)
    .make_tuple_from_list(as.list(x))

as.tuple.set <-
as.tuple.gset <-
as.tuple.cset <-
as.tuple.list <-
function(x)
    do.call(tuple, x)

as.tuple.data.frame <-
function(x)
{
    ret <- as.list(x)
    attributes(ret) <- NULL
    names(ret) <- colnames(x)
    .make_tuple_from_list(ret)
}

as.list.tuple <-
function(x, ...)
    unclass(x)

### cset converters

as.cset <-
function(x)
    UseMethod("as.cset")

as.cset.default <-
function(x)
    cset(gset(x))

as.cset.cset <-
    identity

as.cset.matrix <-
as.cset.data.frame <-
as.cset.logical <-
function(x)
    as.gset(x)

as.cset.tuple <-
as.cset.numeric <-
as.cset.factor <-
as.cset.character <-
as.cset.integer <-
as.cset.list <-
function(x)
{
    sup <- as.gset(x)
    if(!any(duplicated(x))) {
        o <- .list_order(x)
        i <- seq_along(x)
        i[o] <- i
        cset(sup, orderfun = i)
    } else sup
}

as.cset.ordered <-
function(x)
{
    s <- as.character(x)
    o <- .list_order(s)
    dup <- duplicated(s[o])
    cset(as.gset(x), orderfun = order(x)[o][!dup])
}

as.list.cset <-
function(x, ...)
{
    FUN <- cset_orderfun(x)
    L <- .as.list(x)
    ms <- .get_memberships(L)
    if(is.function(FUN) || is.character(FUN)) {
        order <- do.call(FUN, list(L))
        L <- L[order]
        ms <- ms[order]
    } else if(is.integer(FUN) && (length(L) == length(FUN))) {
        L <- L[FUN]
        ms <- ms[FUN]
    }
    .structure(L, memberships = ms)
}

as.character.cset <-
function(x, ...)
{
    x <- as.list(x)
    fac <- unlist(lapply(x, is.factor))
    x[fac] <- lapply(x[fac], as.character)

    FUN <-
        function(x) paste(paste(list(args(x))),
                          paste(as.character(body(x)),
                                collapse = "\n"), sep = "\n")

    fun <- unlist(lapply(x, is.function))
    x[fun] <- lapply(x[fun], FUN)

    paste(x)
}

### make sure that as.list always works

as.list.function <-
function(x, ...)
    list(x)

### .as.list
### In the current implementation, for (g)gets, it's just unclass ...

.as.list <-
function(x, ...)
    UseMethod(".as.list")

.as.list.default <-
function(x, ...)
    as.list(x, ...)

.as.list.set <-
.as.list.gset <-
    unclass

.as.list.cset <-
function(x, ...)
{
    attr(x, "class") <- NULL
    attr(x, "orderfun") <- NULL
    attr(x, "matchfun") <- NULL
    x
}
