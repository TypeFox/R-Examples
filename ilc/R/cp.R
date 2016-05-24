cp <-
function (..., to, as = names, list = character(0), from = 1, 
    o = F) 
{
    names <- sapply(match.call(expand.dots = FALSE)$..., as.character)
    if (length(names) == 0) 
        names <- character(0)
    list <- .Primitive("c")(list, names)
    names <- list
    copied <- rep(T, length(list))
    if (length(as) != length(list)) 
        stop("Wrong number of new names!")
    for (i in 1:length(list)) {
        if (!o && any(find(as[i], numeric = TRUE) %in% to)) {
            copied[i] <- F
        }
        else assign(as[i], get(list[i], pos = from), pos = to)
    }
    if (!sum(copied)) 
        cat("No objects have been copied!\n")
    else cat("Object(s) ", paste(names[copied], ", ", sep = ""), 
        "\n\tcopied from \"", search()[from], "\" to \"", search()[to], 
        "\" as\n", paste(as[copied], ", ", sep = ""), fill = T, 
        sep = "")
    if (!all(copied)) 
        cat("Object(s) ", paste(as[!copied], ", ", sep = ""), 
            "could not be overwritten!\n", fill = T, sep = "")
    invisible(copied)
}
