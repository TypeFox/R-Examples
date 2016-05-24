##' A class for the identifiers of DateTime and ID records in spatial data.
##'
##' The main use of this class and creator function is for
##' \code{\link[sp]{SpatialPointsDataFrame}}s which are used with
##' TimeOrderedRecords for the class \code{trip}.
##'
##' @name TimeOrderedRecords-class
##' @rdname TimeOrderedRecords-class
##' @section Slots:
##'  \describe{
##'    \item{\code{TOR.columns}:}{2-element vector of class \code{"character"}}
##'  }
##' @import sp
##' @import methods
##' @exportClass TimeOrderedRecords
##' @docType class
##' @note  Future versions may change significantly, this class is
##' very basic and could probably be implemented in a better
##' way. Specifying TOR columns by formula would be a useful addition.
##' @seealso \code{\link{TimeOrderedRecords}}, \code{\link{trip}} for creating trip objects, and \code{\link{trip-class}} for that class
##' @examples
##' showClass("TimeOrderedRecords")
##' tor <- new("TimeOrderedRecords", TOR.columns=c("datetime", "ID"))
setClass("TimeOrderedRecords", representation(TOR.columns="character"))

setValidity("TimeOrderedRecords", function(object) {
    if (! is.character(object@TOR.columns) |
        !is.vector(object@TOR.columns)) {
        stop("TimeOrderedRecords data names must be character vector")
        ## also support length == 1?
        if (length(object@TOR.columns) > 2)
            stop("TimeOrderedRecords data names must be of length 2")
        TRUE
    }
})


#' A class for sets of animal trips (track data).
#'
#'
#' An extension of \code{\link[sp]{SpatialPointsDataFrame}} by including
#' \code{"TimeOrderedRecords"}.  The records within the data frame are
#' explicitly ordered by DateTime data within IDs.
#'
#'
#' @name trip-class
#' @rdname trip-class
#' @aliases trip-class lines,trip-method plot,trip,missing-method show,trip-method show,summary.TORdata-method summary,trip-method subset,trip-method
#' @exportClass trip
#' @docType class
#' @section Objects from the Class:
#'
#' Objects can be created by calls of the form
#' \code{trip(obj="SpatialPointsDataFrame", TORnames="TimeOrderedRecords")}.
#' The object contains all the slots present within a
#' \code{\link[sp]{SpatialPointsDataFrame}}, particularly \code{data} which
#' contains columns of at least those specified by \code{TOR.columns}.
#' @seealso
#'
#' \code{\link{trip}} for examples of directly using the class.
#'
#' \code{\link{trip-accessors}} describes methods for accessing information on
#' \code{trip} objects.
#' @keywords classes
#' @examples
#'
#'
#' showClass("trip")
#'
#' ## Examples of general methods
#' ## Continuing the example from '?trip-methods:
#' utils::example("trip-methods", package="trip",
#'                ask=FALSE, echo=FALSE)
#' summary(tr)
#' plot(tr)
#' lines(tr)
#'
#' dim(tr)
#' names(tr)
#' subset(tr, id == "2")
#' as.data.frame(tr)
#'
#' tr[1:3, ]
#' tr[, 1]
#' tr[[1]]
#'
#' if (exists("porpoise")) {
#'     dim(porpoise)
#'     names(porpoise)
#'     porpoise[porpoise[["id"]] == "GUS", ]
#' }
#'
#'
setClass("trip",
         contains=c("TimeOrderedRecords", "SpatialPointsDataFrame"))


##'
#' Internal trip Functions
#'
#' Internal trip functions
#'
#'
#' These are not to be called by the user (or in some cases are just waiting
#' for proper documentation to be written).
#' @name trip-internal
#' @aliases trip-internal .validTORdata .oc.col .single.trip.split .g2ow
#' .gcdist.c .intpFun .tripRbind .distances .abdali
#' @keywords internal
.validTORdata <- function(object) {
    if (!is(object@data, "data.frame"))
        stop("only data frames supported for data slots")
    tid <- as.data.frame(object@data[, object@TOR.columns])
    if (length(tid) == 0)
        stop("timeIDs cannot have zero length")
    if (nrow(tid) < 1)
        stop("no timeIDs set: too few rows")
    if (ncol(tid) < 2)
        stop("no timeIDs set: too few columns")
    if (any(duplicated(as.data.frame(object))))
        stop("duplicated records within data")
    time <- tid[, 1]
    id <- tid[, 2]
    TORlevs <- levels(factor(id))
    if (!is(time, "POSIXt"))
        stop("trip only handles dates and times as POSIXt objects")
    ## mdsumner ID could not be character, because of finite test 2010-04-28
    bad1 <- c(is.na(time), !is.finite(time))
    if (any(bad1))
        return("time data contains missing or non finite values")
    if (any(is.na(id)))
        return("id data contains missing values")
    if (is.numeric(id) & any(!is.finite(id)))
        return("id data contains non-finite values")
    d <- unlist(tapply(time, id, diff))
    if (any(d < 0))
        return("date-times not in order within id")
    if (any(d == 0))
        return("date-times contain duplicates within id")
    short <- which(unlist(tapply(time, id, length)) < 3)
    ## maybe trip enforces this
    if (length(short) > 0) {
        mess <- "\n  fewer than 3 locations for ids:\n"
        mess <- paste(mess,
                      paste(TORlevs[short], collapse=","),
                      sep="")
        return(mess)
    }
    return(TRUE)
}

setValidity("trip", .validTORdata)

## We don't need an S4 class for this, but we do want S4 methods
setOldClass("summary.TORdata")

if (!isClass("ltraj")) setClass("ltraj")
