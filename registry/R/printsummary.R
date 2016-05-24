print.registry <-
function(x, ...)
{
    l <- x$n_of_entries()
    if (l < 1)
        writeLines(gettext("An object of class 'registry' with no entry."))
    else if (l == 1)
        writeLines(gettext("An object of class 'registry' with one entry."))
    else
        writeLines(gettextf("An object of class 'registry' with %d entries.", l))
    invisible(x)
}

print.registry_field <-
function(x, ...) ## remove index function
{
    writeLines(formatUL(x[-7], label = names(x[-7]), ...))
    invisible(x)
}

print.registry_entry <-
function(x, ...)
{
    x <- lapply(.labels(x), paste, collapse = ", ")
    writeLines(formatUL(x, label = names(x)))
    invisible(x)
}

summary.registry <-
function(object, ...)
    as.data.frame(object, ...)

