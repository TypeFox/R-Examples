mv <-
function (..., to, as, list = character(0), from = 1, o = F) 
{
    names <- sapply(match.call(expand.dots = FALSE)$..., as.character)
    if (length(names) == 0) 
        names <- character(0)
    list <- .Primitive("c")(list, names)
    names <- list
    if (missing(as)) 
        as <- names
    copied <- cp(as = as, list = list, to = to, from = from, 
        o = o)
    remove(list = names[copied], pos = from)
    if (sum(copied)) 
        cat("and removed from", from, "\n")
    invisible(copied)
}
