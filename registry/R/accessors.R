"[[.registry" <-
function(x, ...)
{
    if (missing(..1))
        x$get_entry()
    else
        x$get_entry(...)
}

"[.registry" <-
function(x, ...)
{
    if (missing(..1))
        x$get_entries()
    else
        x$get_entries(...)
}

length.registry <-
function(x)
    x$n_of_entries()

