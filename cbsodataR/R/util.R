# Borrow functions from haven, automatically created by 'import_haven_utils.R'
labelled <-
function (x, labels) 
{
    if (!is.numeric(x) && !is.character(x)) {
        stop("`x` must be either numeric or a character vector", 
            call. = FALSE)
    }
    if (typeof(x) != typeof(labels)) {
        stop("`x` and `labels` must be same type", call. = FALSE)
    }
    if (is.null(labels)) {
        stop("`labels` must be a named vector", call. = FALSE)
    }
    structure(x, labels = labels, class = "labelled")
}
as_factor <-
function (x, ...) 
{
    UseMethod("as_factor")
}
as_factor.labelled <-
function (x, levels = c("labels", "values"), ordered = FALSE, 
    ...) 
{
    levels <- match.arg(levels)
    if (is.character(x)) {
        levs <- unname(attr(x, "labels"))
        labs <- switch(levels, labels = names(attr(x, "labels")), 
            values = levs)
        factor(x, levs, labels = labs, ordered = ordered)
    }
    else {
        labs <- attr(x, "labels")
        factor(match(x, labs), levels = unname(labs), labels = names(labs))
    }
}
as_factor.character <-
function (x, ...) 
{
    factor(x, ...)
}
`[.labelled` <-
function (x, ...) 
{
    labelled(NextMethod(), attr(x, "labels"))
}
