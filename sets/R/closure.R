closure <-
function(x, operation, ...)
    UseMethod("closure")

closure.default <-
function(x, operation, ...)
    stop("Not implemented.")

closure.set <-
function(x, operation = c("union", "intersection"), ...)
{
    if (length(x) < 2L) return(x)
    if (!all(sapply(x, is.cset)))
        stop("closure only defined on set of (c,g)sets.")

    if (all(sapply(x, is.set))) {

        dom <- .as.list(do.call(set_union, x))
        x <- lapply(x, .make_list_elements)

        members <-
            binary_closure(do.call(rbind, lapply(x, function(i) dom %in% i)),
                           operation)

        .make_set_from_list(.list_sort(apply(members, 1L,
                                             function(i) .make_set_from_list(dom[i])
                                             )
                                       )
                            )
    } else if (all(sapply(x, is.gset))) {
        operation <- paste("gset_", match.arg(operation), sep = "")
        len <- 0
        while ((newlen <- length(gset_support(x))) != len) {
            len <- newlen
            x <- c(x,
                   as.set(lapply(gset_combn(x, 2L),
                                 function(i) do.call(gset_union, i))
                          )
                   )
        }
        x
    } else {
        operation <- paste("cset_", match.arg(operation), sep = "")
        len <- 0
        while ((newlen <- length(cset_support(x))) != len) {
            len <- newlen
            x <- c(x,
                   as.set(lapply(cset_combn(x, 2L),
                                 function(i) do.call(cset_union, i))
                          )
                   )
        }
        x
    }
}

binary_closure <-
function(x, operation = c("union", "intersection"))
    .Call(R_closure, x,
          pmatch(match.arg(operation), c("union", "intersection")))
