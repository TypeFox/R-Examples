init.CrunchDataset <- function (.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    .Object@variables <- getDatasetVariables(.Object)
    return(.Object)
}
setMethod("initialize", "CrunchDataset", init.CrunchDataset)

getDatasetVariables <- function (x) {
    varcat_url <- variableCatalogURL(x)
    ## Add query params
    return(VariableCatalog(crGET(varcat_url,
        query=list(nosubvars=1, relative="on"))))
}

getNrow <- function (dataset, filtered=TRUE) {
    ## Filtered is the UI "applied_filter". We don't usually want that in R.
    ## TODO: change the default to FALSE
    which.count <- ifelse(isTRUE(filtered), "filtered", "total")

    u <- summaryURL(dataset)
    f <- zcl(activeFilter(dataset))
    q <- crGET(u, query=list(filter=toJSON(f)))
    nrows <- as.integer(round(q$unweighted[[which.count]]))
    return(nrows)
}

##' Is it?
##' @rdname crunch-is
##' @param x an object
##' @return logical
##' @export
is.dataset <- function (x) inherits(x, "CrunchDataset")

##' Name, alias, and description for Crunch objects
##'
##' @param x a Dataset or Variable.
##' @param object Same as \code{x} but for the \code{alias} method, in order to
##' match the generic from another package. Note that \code{alias} is only
##' defined for Variables.
##' @param value For the setters, a length-1 character vector to assign
##' @return Getters return the character object in the specified slot; setters
##' return \code{x} duly modified.
##' @name describe
##' @aliases describe name name<- description description<- alias<- startDate startDate<- endDate endDate<-
##' @seealso \code{\link{Categories}} \code{\link{describe-catalog}}
NULL

##' @rdname describe
##' @export
setMethod("name", "CrunchDataset", function (x) tuple(x)$name)
##' @rdname describe
##' @export
setMethod("name<-", c("CrunchDataset", "character"), function (x, value) {
    out <- setTupleSlot(x, "name", value)
    updateDatasetList() ## could just modify rather than refresh
    invisible(out)
})
##' @rdname describe
##' @export
setMethod("description", "CrunchDataset", function (x) tuple(x)$description)
##' @rdname describe
##' @export
setMethod("description<-", "CrunchDataset", function (x, value) {
    setTupleSlot(x, "description", value)
})
##' @rdname describe
##' @export
setMethod("startDate", "CrunchDataset",
    function (x) trimISODate(tuple(x)$start_date))
##' @rdname describe
##' @export
setMethod("startDate<-", "CrunchDataset", function (x, value) {
    setTupleSlot(x, "start_date", value)
})
##' @rdname describe
##' @export
setMethod("endDate", "CrunchDataset",
    function (x) trimISODate(tuple(x)$end_date))
##' @rdname describe
##' @export
setMethod("endDate<-", "CrunchDataset", function (x, value) {
    setTupleSlot(x, "end_date", value)
})
##' @rdname describe
##' @export
setMethod("id", "CrunchDataset", function (x) tuple(x)$id)

trimISODate <- function (x) {
    ## Drop time from datestring if it's only a date
    if (is.character(x) &&
        nchar(x) > 10 &&
        substr(x, 11, nchar(x)) == "T00:00:00+00:00") {

        x <- substr(x, 1, 10)
    }
    return(x)
}

as.dataset <- function (x, tuple=DatasetTuple()) {
    out <- CrunchDataset(x)
    tuple(out) <- tuple
    return(out)
}

##' Dataset dimensions
##'
##' @param x a Dataset
##' @return integer vector of length 2, indicating the number of rows and
##' non-hidden variables in the dataset. Array subvariables are excluded from
##' the column count.
##' @seealso \code{\link[base]{dim}}
##' @name dim-dataset
NULL

##' @rdname dim-dataset
##' @export
setMethod("dim", "CrunchDataset",
    function (x) c(getNrow(x), ncol(x)))

##' @rdname dim-dataset
##' @export
setMethod("ncol", "CrunchDataset", function (x) length(variables(x)))

namekey <- function (x) {
    if (is.variable(x)) {
        return(match.arg(getOption("crunch.namekey.array"), c("alias", "name")))
    } else {
        return(match.arg(getOption("crunch.namekey.dataset"), c("alias", "name")))
    }
}

##' @rdname describe-catalog
##' @export
setMethod("names", "CrunchDataset", function (x) {
    findVariables(x, key=namekey(x), value=TRUE)
})

##' Dataset weights
##' @param x a Dataset
##' @param value a Variable to set as weight, or NULL to remove the existing
##' weight
##' @return For the getter, a Variable if there is a weight, else NULL. For the
##' setter, x, modified accordingly
##' @export
weight <- function (x) {
    stopifnot(is.dataset(x))
    w <- x@body$weight
    if (!is.null(w)) {
        w <- CrunchVariable(allVariables(x)[[w]], filter=activeFilter(x))
    }
    return(w)
}

##' @rdname weight
##' @export
`weight<-` <- function (x, value) {
    stopifnot(is.dataset(x))
    if (is.variable(value)) {
        value <- self(value)
    } else if (!is.null(value)) {
        halt("Weight must be a Variable or NULL")
    }
    x <- setEntitySlot(x, "weight", value)
    return(x)
}

setMethod("tuple", "CrunchDataset", function (x) x@tuple)
setMethod("tuple<-", "CrunchDataset", function (x, value) {
    x@tuple <- value
    return(x)
})

##' Get a fresh copy from the server
##'
##' Crunch objects usually keep themselves in sync with the server when you
##' manipulate them, but sometimes they can drift. Maybe someone else has
##' modified the dataset you're working on, or maybe
##' you have modified a variable outside of the context of its dataset.
##' refresh() allows you to get back in sync.
##'
##' @param x pretty much any Crunch object
##' @return a new version of \code{x}
##' @name refresh
##' @aliases refresh
NULL

##' @rdname refresh
##' @export
setMethod("refresh", "CrunchDataset", function (x) {
    tup <- refresh(tuple(x))
    out <- as.dataset(crGET(self(x)), tuple=tup)
    duplicates(allVariables(out)) <- duplicates(allVariables(x))
    activeFilter(out) <- activeFilter(x)
    return(out)
})

##' Delete a Crunch object from the server
##'
##' These methods delete entities, notably Datasets and Variables within them,
##' from the server. This action is permanent and cannot be undone, so it
##' should not be done lightly. Consider instead using \code{archive}
##' for datasets and \code{\link{hide}} for variables
##'
##' @param x a Crunch object
##' @param confirm logical: should the user be asked to confirm deletion.
##' Option available for datasets and teams only. Default is \code{TRUE} if in
##' an interactive session. You can avoid the confirmation prompt if you delete
##' \code{with(\link{consent})}.
##' @param ... additional arguments, in the generic
##' @seealso \code{\link{hide}} \code{\link{deleteDataset}}
##' @name delete
##' @aliases delete
NULL

##' @rdname delete
##' @export
setMethod("delete", "CrunchDataset",
    function (x, confirm=requireConsent() | is.readonly(x), ...) {
        out <- delete(tuple(x), confirm=confirm)
        invisible(out)
    })

##' @export
as.list.CrunchDataset <- function (x, ...) {
    lapply(seq_along(variables(x)), function (i) x[[i]])
}

##' See the appended batches of this dataset
##' @param x a \code{CrunchDataset}
##' @return a \code{BatchCatalog}
##' @export
batches <- function (x) BatchCatalog(crGET(shojiURL(x, "catalogs", "batches")))

joins <- function (x) ShojiCatalog(crGET(shojiURL(x, "catalogs", "joins")))

setDatasetVariables <- function (x, value) {
    v <- urls(value)
    x@variables[v] <- value
    ordering(x@variables) <- ordering(value)
    return(x)
}

setMethod("datasetReference", "CrunchDataset", function (x) self(x))

variableCatalogURL <- function (dataset) {
    ## Get the variable catalog URL that corresponds to an object
    if (class(dataset) == "VariableCatalog") return(self(dataset))
    if (!is.dataset(dataset)) {
        dataset <- ShojiObject(crGET(datasetReference(dataset)))
    }
    return(shojiURL(dataset, "catalogs", "variables"))
}

summaryURL <- function (x) shojiURL(x, "views", "summary")

##' Access a Dataset's Variables Catalog
##'
##' Datasets contain collections of variables. For a few purposes, such as
##' editing variables' metadata, it is helpful to access these variable catalogs
##' more directly.
##'
##' \code{variables} gives just the active variables in the dataset, while
##' \code{allVariables}, as the name suggests, yields all variables, including
##' hidden variables.
##' @param x a Dataset
##' @param value For the setters, a VariableCatalog to assign.
##' @return Getters return VariableCatalog; setters return \code{x} duly
##' modified.
##' @name dataset-variables
##' @aliases dataset-variables variables variables<- allVariables allVariables<-
NULL

##' @rdname dataset-variables
##' @export
setMethod("variables", "CrunchDataset", function (x) active(allVariables(x)))
##' @rdname dataset-variables
##' @export
setMethod("variables<-", c("CrunchDataset", "VariableCatalog"),
    setDatasetVariables)
##' @rdname dataset-variables
##' @export
setMethod("allVariables", "CrunchDataset", function (x) x@variables)
##' @rdname dataset-variables
##' @export
setMethod("allVariables<-", c("CrunchDataset", "VariableCatalog"),
    setDatasetVariables)

setMethod("hidden", "CrunchDataset", function (x) hidden(allVariables(x)))


webURL <- function (x) {
    ## URL to view this dataset in the web app
    stopifnot(is.dataset(x))
    return(paste0(absoluteURL("/", getOption("crunch.api")), "dataset/", id(x)))
}

##' as.environment method for CrunchDataset
##'
##' This method allows you to \code{eval} within a Dataset.
##'
##' @param x CrunchDataset
##' @return an environment in which named objects are (promises that return)
##' CrunchVariables.
setMethod("as.environment", "CrunchDataset", function (x) {
    out <- new.env()
    out$.crunchDataset <- x
    with(out, for (v in aliases(allVariables(x))) eval(parse(text=paste0("delayedAssign('", v, "', .crunchDataset[['", v, "']])"))))
    return(out)
})

.releaseDataset <- function (dataset) {
    release_url <- absoluteURL("release/", self(dataset))
    crPOST(release_url)
}
