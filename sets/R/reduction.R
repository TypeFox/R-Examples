reduction <-
function(x, operation, ...)
    UseMethod("reduction")

reduction.default <-
function(x, operation, ...)
    stop("Not implemented.")

reduction.set <-
function(x, operation = c("union", "intersection"), ...)
{
    operation <- match.arg(operation)
    if (length(x) < 2L) return(x)
    if (!all(sapply(x, is.cset)))
        stop("reduction only defined on set of (c,g)sets.")

    if (all(sapply(x, is.set))) {

        dom <- .as.list(do.call(set_union, x))
        x <- lapply(x, .make_list_elements)

        members <-
            binary_reduction(do.call(rbind, lapply(x, function(i) dom %in% i)),
                             operation)

        .make_set_from_list(.list_sort(apply(members, 1L,
                                             function(i) .make_set_from_list(dom[i])
                                             )
                                       )
                            )
    } else if (all(sapply(x, is.gset))) {
        clo <- closure(x, operation)
        FUN <- function(e)
            !isTRUE(gset_is_equal(closure(gset_difference(x, set(e)),
                                          operation), clo))
        as.set(Filter(FUN, .as.list(x)))
    } else {
        clo <- closure(x, operation)
        FUN <- function(e)
            !isTRUE(cset_is_equal(closure(cset_difference(x, set(e)),
                                          operation), clo))
        as.set(Filter(FUN, .as.list(x)))
    }
}

binary_reduction <-
function(x, operation = c("union", "intersection"))
    .Call(R_reduction, x,
          pmatch(match.arg(operation), c("union", "intersection")))
