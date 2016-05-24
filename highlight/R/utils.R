# because utils:::SweaveParseOptions is now forbidden on CRAN
SweaveParseOptions <- function (text, defaults = list(), check = NULL){
    x <- sub("^[[:space:]]*(.*)", "\\1", text)
    x <- sub("(.*[^[:space:]])[[:space:]]*$", "\\1", x)
    x <- unlist(strsplit(x, "[[:space:]]*,[[:space:]]*"))
    x <- strsplit(x, "[[:space:]]*=[[:space:]]*")
    if (length(x)) {
        if (length(x[[1L]]) == 1L)
            x[[1L]] <- c("label", x[[1L]])
    }
    else return(defaults)
    if (any(sapply(x, length) != 2L))
        stop(gettextf("parse error or empty option in\n%s", text),
            domain = NA)
    options <- defaults
    for (k in seq_along(x)) options[[x[[k]][1L]]] <- x[[k]][2L]
    if (!is.null(options[["label"]]) && !is.null(options[["engine"]]))
        options[["label"]] <- sub(paste0("\\.", options[["engine"]],
            "$"), "", options[["label"]])
    if (!is.null(check))
        check(options)
    else options
}
