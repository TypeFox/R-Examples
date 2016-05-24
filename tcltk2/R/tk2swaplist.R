tk2swaplist <- function(items, selection, title = "Select items", ...)
{
    win <- tktoplevel()
    res <- try(tclRequire("swaplist"), silent = TRUE)
    if (inherits(res, "try-error"))
        stop("swaplist Tcl package not available")
    sel <- tclVar()
    res <- tcl("swaplist::swaplist", win, sel, items, selection,
        title = title, ...)
    if (tclvalue(res) == 0) { # User cancelled
        res <- character(0)
    } else res <- tclObj(sel)
    if (is.ordered(items))
        return(ordered(as.character(res), levels = levels(items)))
    if (is.factor(items))
        return(factor(as.character(res), levels = levels(items)))
    switch(typeof(items),
        integer = as.integer(res),
        double = as.numeric(res),
        logical = as.logical(res),
        complex = as.complex(res),
        as.character(res))
}
