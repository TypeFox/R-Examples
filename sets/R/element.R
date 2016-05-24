## * Element function, used for the creation of gsets by elements

e <-
function(x, memberships = 1L)
{
    if (is_element(x)) return(x)
    if (length.set(memberships) > 1L)
        memberships <- as.gset(memberships)
    else if (isTRUE(memberships >= 1))
        memberships <- as.integer(memberships)
    .stop_if_memberships_are_invalid(memberships)
    .make_element_from_support_and_memberships(x, memberships)
}

.make_element_from_support_and_memberships <-
function(x, memberships)
    .structure(list(x),
               memberships = memberships,
               class = "element")


print.element <-
function(x, ...)
{
    writeLines(format(x, ...))
    invisible(x[[1]])
}

format.element <-
function(x, ...)
{
    paste(paste(LABEL(x[[1]], ...), collapse = " "),
          " [",
          paste(format(.get_memberships(x)), collapse = ", "),
          "]",
          sep = "")
}

LABEL.element <-
function(x, limit, ...)
    format(x, ...)

