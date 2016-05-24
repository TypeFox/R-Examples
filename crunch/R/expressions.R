##' Construct Crunch Expressions
##'
##' Crunch Expressions, i.e. \code{CrunchExpr} and \code{CrunchLogicalExpr},
##' encapuslate derivations of Crunch variables, which are only evaluated when
##' passed to a function like \code{as.vector}. They allow you to compose
##' functional expressions of variables and evaluate them against the server
##' only when appropriate.
##'
##' @param x an input
##' @param table For \code{\%in\%}. See \code{\link[base]{match}}
##' @param resolution For \code{rollup}. Either \code{NULL} or a character in
##' c("Y", "Q", "M", "W", "D", "h", "m", "s", "ms") indicating the unit of
##' time at which a Datetime variable should be aggregated. If \code{NULL},
##' the server will determine an appropriate resolution based on the range of
##' the data.
##' @return Most functions return a CrunchExpr or CrunchLogicalExpr.
##' \code{as.vector} returns an R vector.
##' @aliases expressions
##' @name expressions
NULL

##' @rdname variable-to-R
##' @export
setMethod("as.vector", "CrunchExpr", function (x, mode) {
    payload <- list(command="select", variables=list(out=zcl(x)))
    if (length(x@filter)) {
        payload[["filter"]] <- x@filter
    }
    out <- paginatedGET(paste0(x@dataset_url, "table/"),
        query=list(query=toJSON(payload)),
        limit=10000000000) ## Because server doesn't accept pagination yet
    ## pass in the variable metadata to the column parser
    variable <- VariableEntity(structure(list(body=out$metadata$out),
        class="shoji"))
    return(columnParser(out$metadata$out$type, mode)(out$data$out, variable))
})

##' @rdname toVariable
##' @export
setMethod("toVariable", "CrunchExpr", function (x, ...) {
    structure(list(expr=zcl(x), ...), class="VariableDefinition")
})

## "Ops" for Crunch Variables
##
## Most of the indirection here is to programatically create the Ops methods
## for the right combinations of multiple-dispatch signatures

math.exp <- function (e1, e2, operator) {
    ## Generic function that creates CrunchExpr of `e1 %operator% e2`
    ex <- zfunc(operator, e1, e2)
    ds.url <- unique(unlist(lapply(list(e1, e2), datasetReference))) %||% ""
    logics <- c("in", "<", ">", ">=", "<=", "==", "!=", "and", "or")
    if (operator %in% logics) {
        Constructor <- CrunchLogicalExpr
    } else {
        Constructor <- CrunchExpr
    }
    ## TODO: if e1 and e2 have filters, pass them along, and if both do,
    ## make sure that they're the same
    return(Constructor(expression=ex, dataset_url=ds.url))
}

vxr <- function (i) {
    ## Create math.exp of Variable x R.object
    force(i)
    return(function (e1, e2) math.exp(e1, typeof(e2, e1), i))
}

rxv <- function (i) {
    ## Create math.exp of R.object x Variable
    force(i)
    return(function (e1, e2) math.exp(typeof(e1, e2), e2, i))
}

vxv <- function (i) {
    ## Create math.exp of two non-R.objects
    force(i)
    return(function (e1, e2) math.exp(e1, e2, i))
}

.sigs <- list(
    c("TextVariable", "character"),
    c("NumericVariable", "numeric"),
    c("DatetimeVariable", "Date"),
    c("DatetimeVariable", "POSIXt"),
    c("CategoricalVariable", "numeric")
)

.rtypes <- unique(vapply(.sigs, function (a) a[[2]], character(1)))
.nomath <- which(!vapply(.sigs,
    function (a) a[[1]] %in% c("TextVariable", "CategoricalVariable"),
    logical(1)))

for (i in c("+", "-", "*", "/", "<", ">", ">=", "<=")) {
    for (j in .nomath) {
        setMethod(i, .sigs[[j]], vxr(i))
        setMethod(i, rev(.sigs[[j]]), rxv(i))
    }
    for (j in setdiff(.rtypes, "character")) {
        setMethod(i, c("CrunchExpr", j), vxv(i)) # no typeof?
        setMethod(i, c(j, "CrunchExpr"), vxv(i)) # no typeof?
    }
    setMethod(i, c("CrunchVariable", "CrunchVariable"), vxv(i))
    setMethod(i, c("CrunchExpr", "CrunchVariable"), vxv(i))
    setMethod(i, c("CrunchVariable", "CrunchExpr"), vxv(i))
}

.catmeth <- function (i, Rarg=1) {
    force(i)
    force(Rarg)
    return(function (e1, e2) {
        if (Rarg == 1) {
            e1 <- typeof(n2i(as.character(e1), categories(e2)), e2)
        } else {
            e2 <- typeof(n2i(as.character(e2), categories(e1)), e1)
        }
        return(math.exp(e1, e2, i))
    })
}

for (i in c("==", "!=")) {
    for (j in seq_along(.sigs)) {
        setMethod(i, .sigs[[j]], vxr(i))
        setMethod(i, rev(.sigs[[j]]), rxv(i))
    }
    setMethod(i, c("CategoricalVariable", "character"), .catmeth(i, 2))
    setMethod(i, c("CategoricalVariable", "factor"), .catmeth(i, 2))
    setMethod(i, c("character", "CategoricalVariable"), .catmeth(i, 1))
    setMethod(i, c("factor", "CategoricalVariable"), .catmeth(i, i))
    for (j in .rtypes) {
        setMethod(i, c("CrunchExpr", j), vxv(i)) # no typeof?
        setMethod(i, c(j, "CrunchExpr"), vxv(i)) # no typeof?
    }
    setMethod(i, c("CrunchVariable", "CrunchVariable"), vxv(i))
    setMethod(i, c("CrunchExpr", "CrunchVariable"), vxv(i))
    setMethod(i, c("CrunchVariable", "CrunchExpr"), vxv(i))
}

setMethod("&", c("CrunchExpr", "CrunchExpr"), vxv("and"))
setMethod("&", c("logical", "CrunchExpr"), vxv("and"))
setMethod("&", c("CrunchExpr", "logical"), vxv("and"))
setMethod("|", c("CrunchExpr", "CrunchExpr"), vxv("or"))
setMethod("|", c("logical", "CrunchExpr"), vxv("or"))
setMethod("|", c("CrunchExpr", "logical"), vxv("or"))

##' @rdname expressions
##' @export
setMethod("!", c("CrunchExpr"),
    function (x) {
        CrunchLogicalExpr(expression=zfunc("not", x),
            dataset_url=datasetReference(x))
    })

##' @importFrom utils head tail
.seqCrunch <- function (x, table) {
    ## Given x %in% table, if table is numeric, see if we can/should collapse
    ## it into a range query rather than sending lots of distinct values

    if (is.numeric(table) &&
        length(table) > 2 &&
        identical(as.numeric(head(table, 1):tail(table, 1)), as.numeric(table))) {

        return(zfunc("between",
            x,
            head(table, 1),
            tail(table, 1) + 1))
            ## Add 1 because "between" by default doesn't include the upper
            ## bound and explicitly overriding that is failing. See #112089103.
            ## When that is fixed, we can do the following:
            #rep(TRUE, 2L))) ## Inclusive on both sides
    } else {
        return(zfunc(ifelse(length(table) == 1L, "==", "in"), x, table))
    }
}

.inCrunch <- function (x, table) {
    ## TODO: bring in .seqCrunch here, once it is working
    math.exp(x, typeof(table, x), ifelse(length(table) == 1L, "==", "in"))
}

##' @rdname expressions
##' @export
setMethod("%in%", c("CategoricalVariable", "character"),
    function (x, table) .inCrunch(x, n2i(table, categories(x))))
##' @rdname expressions
##' @export
setMethod("%in%", c("CategoricalVariable", "factor"),
    function (x, table) x %in% as.character(table))

## Iterated version of below:
for (i in seq_along(.sigs)) {
    setMethod("%in%", .sigs[[i]], .inCrunch)
}

##' @rdname expressions
##' @export
setMethod("%in%", c("TextVariable", "character"), .inCrunch)
##' @rdname expressions
##' @export
setMethod("%in%", c("NumericVariable", "numeric"), .inCrunch)
##' @rdname expressions
##' @export
setMethod("%in%", c("DatetimeVariable", "Date"), .inCrunch)
##' @rdname expressions
##' @export
setMethod("%in%", c("DatetimeVariable", "POSIXt"), .inCrunch)
##' @rdname expressions
##' @export
setMethod("%in%", c("CategoricalVariable", "numeric"), .inCrunch)

setMethod("datasetReference", "CrunchExpr", function (x) x@dataset_url)

##' @rdname expressions
##' @export
setMethod("is.na", "CrunchVariable", function (x) {
    CrunchLogicalExpr(expression=zfunc("is_missing", x),
        dataset_url=datasetReference(x))
})

##' @rdname expressions
##' @export
bin <- function (x) {
    ## TODO: this should be filtered by x's filter
    CrunchExpr(expression=zfunc("bin", x),
        dataset_url=datasetReference(x) %||% "")
}


##' @rdname expressions
##' @export
rollup <- function (x, resolution=rollupResolution(x)) {
    valid_res <- c("Y", "Q", "M", "W", "D", "h", "m", "s", "ms")
    force(resolution)
    if (!is.null(resolution) && !(resolution %in% valid_res)) {
        halt(dQuote("resolution"), " is invalid. Valid values are NULL, ",
            serialPaste(valid_res))
    }
    if (is.variable(x) && !is.Datetime(x)) {
        halt("Cannot rollup a variable of type ", dQuote(type(x)))
    }
    ## TODO: this should be filtered by x's filter
    CrunchExpr(expression=zfunc("rollup", x, list(value=resolution)),
        ## list() so that the resolution value won't get typed
        dataset_url=datasetReference(x) %||% "")
}

rollupResolution <- function (x) {
    if (is.Datetime(x)) {
        return(tuple(x)$rollup_resolution)
    } else {
        return(NULL)
    }
}

##' @rdname variable-extract
##' @export
setMethod("[", c("CrunchExpr", "CrunchLogicalExpr"), function (x, i, ...) {
    f <- activeFilter(x)
    if (length(zcl(f))) {
        i <- f & i
    }
    activeFilter(x) <- i
    return(x)
})
