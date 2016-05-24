parse_column <- list(
    numeric=function (col, variable) {
        missings <- vapply(col, Negate(is.numeric), logical(1))
        col[missings] <- NA_real_
        return(as.numeric(unlist(col)))
    },
    text=function (col, variable) {
        missings <- vapply(col, Negate(is.character), logical(1))
        col[missings] <- NA_character_
        return(as.character(unlist(col)))
    },
    categorical=function (col, variable) {
        out <- columnParser("numeric")(col)
        cats <- na.omit(categories(variable))
        out <- factor(names(cats)[match(out, ids(cats))], levels=names(cats))
        return(out)
    },
    categorical_ids=function (col, variable) {
        missings <- vapply(col, is.list, logical(1)) ## for the {?:values}
        col[missings] <- lapply(col[missings], function (x) x[["?"]])
        return(as.numeric(unlist(col)))
    },
    categorical_numeric_values=function (col, variable) {
        out <- columnParser("numeric")(col)
        cats <- na.omit(categories(variable))
        out <- values(cats)[match(out, ids(cats))]
        return(out)
    },
    categorical_array=function (col, variable) {
        out <- columnParser("categorical")(unlist(col), variable)
        ncols <- length(tuple(variable)$subvariables)
        out <- t(structure(out, .Dim=c(ncols, length(out)/ncols)))
        return(out)
    },
    datetime=function (col, variable) {
        out <- columnParser("text")(col)
        if (all(grepl("[0-9]{4}-[0-9]{2}-[0-9]{2}", out))) {
            ## return Date if resolution >= D
            return(as.Date(out))
        } else {
            ## TODO: use from8601, defined below
            return(as.POSIXct(out))
        }
    }
)
columnParser <- function (vartype, mode=NULL) {
    if (vartype == "categorical") {
        ## Deal with mode. Valid modes: factor (default), numeric, id
        if (!is.null(mode)) {
            if (mode == "numeric") {
                vartype <- "categorical_numeric_values"
            } else if (mode == "id") {
                vartype <- "categorical_ids" ## The numeric parser will return ids, right?
            }
        }
    }
    return(parse_column[[vartype]] %||% parse_column[["numeric"]])
}

getValues <- function (x, ...) {
    paginatedGET(shojiURL(x, "views", "values"), list(...))
}

paginatedGET <- function (url, query, offset=0,
                          limit=getOption("crunch.page.size") %||% 1000) {
    ## Paginate the GETting of values. Called both from getValues and in
    ## the as.vector.CrunchExpr method in expressions.R

    query$offset <- offset
    query$limit <- limit

    out <- list()
    keep.going <- TRUE
    i <- 1
    with(temp.option(scipen=15), {
        ## Mess with scipen so that the query string formatter doesn't
        ## convert an offset like 100000 to '1+e05', which server rejects
        while(keep.going) {
            out[[i]] <- crGET(url, query=query)
            if (length(out[[i]]) < limit) {
                keep.going <- FALSE
            } else {
                query$offset <- query$offset + limit
                i <- i + 1
            }
        }
    })
    return(unlist(out, recursive=FALSE))
}

##' Convert Variables to local R objects
##'
##' @param x a CrunchVariable subclass
##' @param mode for Categorical variables, one of either "factor" (default,
##' which returns the values as factor); "numeric" (which returns the numeric
##' values); or "id" (which returns the category ids). If "id", values
##' corresponding to missing categories will return as the underlying integer
##' codes; i.e., the R representation will not have any \code{NA}s. Otherwise,
##' missing categories will all be returned \code{NA}. For non-Categorical
##' variables, the \code{mode} argument is ignored.
##' @return an R vector of the type corresponding to the Variable. E.g.
##' CategoricalVariable yields type factor by default, NumericVariable yields
##' numeric, etc.
##' @name variable-to-R
NULL

##' @rdname variable-to-R
##' @export
setMethod("as.vector", "CrunchVariable", function (x, mode) {
    f <- zcl(activeFilter(x))
    columnParser(type(x), mode)(getValues(x, filter=toJSON(f)), x)
})

from8601 <- function (x) {
    ## Crunch timestamps look like "2015-02-12T10:28:05.632000+00:00"

    ## TODO: pull out the ms, as.numeric them, and add to the parsed date
    ## Important for the round trip of datetime data

    ## First, strip out ms and the : in the time zone
    x <- sub("\\.[0-9]+", "", sub("^(.*[+-][0-9]{2}):([0-9]{2})$", "\\1\\2", x))
    ## Then parse
    return(strptime(x, "%Y-%m-%dT%H:%M:%S%z"))
}
