joinDatasets <- function (x, y, by=intersect(names(x), names(y)), by.x=by, 
                         by.y=by, all=FALSE, all.x=TRUE, all.y=FALSE) {

    ## Validate inputs
    if (!is.dataset(x)) halt("x must be a Crunch Dataset")
    if (!is.dataset(y)) halt("y must be a Crunch Dataset")
    if (length(by.x) > 1 || length(by.y) > 1) {
        halt("Can only join 'by' a single key")
    }
    if (is.numeric(by.x)) by.x <- names(x)[by.x]
    if (is.numeric(by.y)) by.y <- names(y)[by.y]

    ## Get var urls for by.x and by.y
    x_url <- findVariableURLs(x, by.x)
    if (length(x_url) == 0) {
        halt("Join key not found in first Dataset")
    }
    if (length(x_url) > 1) {
        halt("Multiple variables matched ", dQuote(by.x),
            " in ", namekey(x), " of first Dataset")
    }
    y_url <- findVariableURLs(y, by.y)
    if (length(y_url) == 0) {
        halt("Join key not found in second Dataset")
    }
    if (length(y_url) > 1) {
        halt("Multiple variables matched ", dQuote(by.y),
            " in ", namekey(y), " of second Dataset")
    }
    ## Get join catalog url
    join_url <- x@catalogs$joins

    payload <- structure(list(list(left_key=x_url, right_key=y_url)),
        .Names=paste0(join_url, tuple(y)$id, "/"))
    crPATCH(join_url, body=toJSON(payload))
    invisible(refresh(x))
}
