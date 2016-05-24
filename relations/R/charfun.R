## Characteristic functions.

### * relation_charfun

relation_charfun <-
function(x, components = FALSE)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    local({I <- .incidence(x)
           D <- .domain(x)
           a <- .arity(x)
           f <- if (!components)
               function(...) {
                   args <- list(...)
                   if (a == 2L) {
                       ## recycle for binary relations
                       maxlen <- max(sapply(args, length))
                       args <- lapply(args, rep, length.out = maxlen)
                   }
                   if(length(args) != a)
                       stop("Wrong number of arguments.")
                   t <- .split_into_components(do.call("cbind", args))
                   ind <- rbind(mapply(.exact_match, t, D))
                   if (any(is.na(ind)))
                       stop("Out-of-domain elements.")
                   as.logical(I[ind])
               }
           else
               function(t) {
                   t <- .split_into_components(t)
                   ind <- rbind(mapply(.exact_match, t, D))
                   if (any(is.na(ind)))
                       stop("Out-of-domain elements.")
                   as.logical(I[ind])
               }
           class(f) <- "relation_charfun"
           f
       })
}

print.relation_charfun <-
function(x, ...)
{
    writeLines(gettextf("The characteristic function of a relation with arity %d.",
                        environment(x)$a))
    invisible(x)
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
