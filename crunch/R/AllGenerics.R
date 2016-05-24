##' @export
setGeneric("values", function (x) standardGeneric("values"))
setGeneric("values<-", function (x, value) standardGeneric("values<-"))

setGeneric("id", function (x) standardGeneric("id"))
setGeneric("is.selected", function (x) standardGeneric("is.selected"))

setGeneric("ids", function (x) standardGeneric("ids"))
setGeneric("ids<-", function (x, value) standardGeneric("ids<-"))
setGeneric("is.dichotomized", function (x) standardGeneric("is.dichotomized"))
setGeneric("dichotomize", function (x, i) standardGeneric("dichotomize"))
setGeneric("undichotomize", function (x) standardGeneric("undichotomize"))
setGeneric("value", function (x) standardGeneric("value"))
setGeneric("value<-", function (x, value) standardGeneric("value<-"))
setGeneric("name", function (x) standardGeneric("name"))
setGeneric("name<-", function (x, value) standardGeneric("name<-"))
setGeneric("expr", function (x) standardGeneric("expr"))
setGeneric("expr<-", function (x, value) standardGeneric("expr<-"),
    signature="x")
setGeneric("description", function (x) standardGeneric("description"))
setGeneric("description<-",
    function (x, value) standardGeneric("description<-"), signature="x")
setGeneric("startDate", function (x) standardGeneric("startDate"))
setGeneric("startDate<-",
function (x, value) standardGeneric("startDate<-"), signature="x")
setGeneric("endDate", function (x) standardGeneric("endDate"))
setGeneric("endDate<-",
    function (x, value) standardGeneric("endDate<-"), signature="x")
setGeneric("alias")
setGeneric("alias<-", function (x, value) standardGeneric("alias<-"),
    signature="x")
setGeneric("aliases", function (x) standardGeneric("aliases"))
setGeneric("aliases<-", function (x, value) standardGeneric("aliases<-"),
    signature="x")
setGeneric("descriptions", function (x) standardGeneric("descriptions"))
setGeneric("descriptions<-",
    function (x, value) standardGeneric("descriptions<-"), signature="x")
setGeneric("emails", function (x) standardGeneric("emails"))
setGeneric("types", function (x) standardGeneric("types"))
setGeneric("timestamps", function (x) standardGeneric("timestamps"))

setGeneric("type", function (x) standardGeneric("type"))
setGeneric("type<-", function (x, value) standardGeneric("type<-"))

setGeneric("categories", function (x) standardGeneric("categories"))
setGeneric("categories<-", function (x, value) standardGeneric("categories<-"))
setGeneric("variables", function (x) standardGeneric("variables"))
setGeneric("variables<-", function (x, value) standardGeneric("variables<-"))
setGeneric("allVariables", function (x) standardGeneric("allVariables"))
setGeneric("allVariables<-",
    function (x, value) standardGeneric("allVariables<-"))
setGeneric("subvariables", function (x) standardGeneric("subvariables"))
setGeneric("subvariables<-",
    function (x, value) standardGeneric("subvariables<-"))
setGeneric("datasetReference", function (x) standardGeneric("datasetReference"))
setGeneric("hide", function (x) standardGeneric("hide"))
setGeneric("unhide", function (x) standardGeneric("unhide"))

setGeneric("urls", function (x) standardGeneric("urls"))
setGeneric("self", function (x) standardGeneric("self"))
setGeneric("refresh", function (x) standardGeneric("refresh"))
setGeneric("delete", function (x, ...) standardGeneric("delete"),
    signature="x")
setGeneric("readonly<-", function (x, value) standardGeneric("readonly<-"))
setGeneric("entities", function (x, ...) standardGeneric("entities"))
setGeneric("entities<-", function (x, value) standardGeneric("entities<-"))
setGeneric("tuple", function (x) standardGeneric("tuple"))
setGeneric("tuple<-", function (x, value) standardGeneric("tuple<-"))
setGeneric("ordering", function (x) standardGeneric("ordering"))
setGeneric("ordering<-", function (x, value) standardGeneric("ordering<-"))
setGeneric("duplicates", function (x) standardGeneric("duplicates"))
setGeneric("duplicates<-", function (x, value) standardGeneric("duplicates<-"))
setGeneric("entity", function (x) standardGeneric("entity"))
setGeneric("index", function (x) standardGeneric("index"))
setGeneric("index<-", function (x, value) standardGeneric("index<-"))
setGeneric("active", function (x) standardGeneric("active"))
setGeneric("hidden", function (x) standardGeneric("hidden"))
setGeneric("archived", function (x) standardGeneric("archived"))
setGeneric("imported", function (x) standardGeneric("imported"))
setGeneric("pending", function (x) standardGeneric("pending"))
setGeneric("permissions", function (x) standardGeneric("permissions"))
setGeneric("members", function (x) standardGeneric("members"))
setGeneric("members<-", function (x, value) standardGeneric("members<-"))
setGeneric("filters", function (x) standardGeneric("filters"))
setGeneric("filters<-", function (x, value) standardGeneric("filters<-"))
setGeneric("appliedFilters", function (x) standardGeneric("appliedFilters"))
setGeneric("appliedFilters<-",
    function (x, value) standardGeneric("appliedFilters<-"))
setGeneric("activeFilter", function (x) standardGeneric("activeFilter"))
setGeneric("activeFilter<-",
    function (x, value) standardGeneric("activeFilter<-"))
setGeneric("is.public", function (x) standardGeneric("is.public"))
setGeneric("is.public<-", function (x, value) standardGeneric("is.public<-"))

setGeneric("dim")
setGeneric("ncol")
setGeneric("mean")
setGeneric("length")
setGeneric("sd")
setGeneric("median")
setGeneric("min")
setGeneric("max")
setGeneric("na.omit")
setGeneric("as.vector")
setGeneric("as.environment")
setGeneric("dimnames")
setGeneric("margin.table")

setGeneric("prop.table")
setGeneric("round")

setGeneric("subset")

##' Generic method for converting objects to Crunch representations
##'
##' If you have other object types you wish to convert to Crunch variables,
##' you can declare methods for \code{toVariable}
##' @param x the object
##' @param ... additional arguments
##' @return a list object suitable for POSTing to the Crunch API. See the API
##' documentation for specifications.
##' @rdname toVariable
##' @aliases toVariable
##' @export
setGeneric("toVariable", function (x, ...) standardGeneric("toVariable"))

setGeneric("lapply")
setGeneric("is.na")
setGeneric("is.na<-")
setGeneric("%in%")

setGeneric("zcl", function (x) standardGeneric("zcl"))

##' toJSON methods for Crunch objects
##'
##' \code{crunch} uses the \code{jsonlite} package for (de)serialization of
##' JSON. Unlike \code{RJSONIO}'s \code{toJSON}, \code{\link[jsonlite]{toJSON}}
##' does not allow for defining S4 methods for other object types. So,
##' \code{crunch::toJSON} wraps \code{jsonprep}, which exists to translate
##' objects to base R objects, which \code{jsonlite::toJSON} can handle.
##' \code{jsonprep} is defined as an S4 generic, and it is exported (unlike
##' code{jsonlite::asJSON}), so you can define methods for it if you have other
##' objects that you want to successfully serialize to JSON.
##'
##' @param x the object
##' @param ... additional arguments
##' @return \code{jsonprep} returns a base R object that \code{jsonlite::toJSON}
##' can handle. \code{toJSON} returns the JSON-serialized character object.
##' @name tojson-crunch
##' @seealso \code{\link[jsonlite]{toJSON}}
NULL

##' @rdname tojson-crunch
##' @export
setGeneric("jsonprep", function (x, ...) standardGeneric("jsonprep"))

setGeneric("getShowContent",
    function (x, ...) standardGeneric("getShowContent"))
