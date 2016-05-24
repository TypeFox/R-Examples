#' @include cfData.R cfStation.R

# cfDataList --------------------------------------------------------------

#' @importFrom methods setClass
setClass("cfDataList", contains = "list")

# Methods -----------------------------------------------------------------

#' @importFrom methods setClass slot
setMethod("show",
          signature(object = "cfDataList"),
          function (object)
          {
            n_row = sapply(object@.Data, nrow)
            data = sapply(object, slot, "dt_name")
            type = sapply(object, slot, "dt_type")
            type[is.na(type)] = ""
            start_list = lapply(object, function(x)
              min(as(x, "data.frame")[, 2]))
            start = sapply(start_list, format, "(%Y-%m-%d %k:00)")
            end_list = lapply(object, function(x)
              max(as(x, "data.frame")[, 2]))
            end = sapply(end_list, format, "(%Y-%m-%d %k:00)")
            cat("List containing clifro data frames:\n")
            print(data.frame(data = tidy_names(data, max_len = 15),
                             type = type,
                             start = start,
                             end = end,
                             rows = n_row,
                             row.names = paste0("df ", seq_along(n_row), ")")))
          }
)

#' Subsetting Methods for Clifro Objects
#'
#' Operators acting on \code{cfDataList}, \code{cfDatatype}, \code{cfStation},
#' and \code{dataFrame} objects.
#'
#' These are methods for the generic operators for classes within \pkg{clifro}.
#' They are intended to give the user the familiar functionality of subsetting
#' \code{\link{data.frame}} objects.
#'
#' @param x a \pkg{clifro} object
#' @param i indices specifying elements to extract. Indices are
#'          \code{numeric} or \code{character} vectors or empty (missing) or
#'          \code{NULL}. Character vectors will be matched to the names of
#'          the object.
#' @param j indices specifying elements to extract. Indices are
#'          \code{numeric} or \code{character} vectors or empty (missing) or
#'          \code{NULL}. Character vectors will be matched to the names of
#'          the object.
#' @param name a literal character string. This is partially matched to the
#'             names of the object.
#' @param drop if \code{TRUE}, the result is coerced to the lowest possible
#'             dimension. See \code{\link{drop}} for further details.
#'
#' @docType methods
#' @name Extract
#' @rdname Extract
#' @aliases [,cfDataList,ANY,ANY,ANY-method
#'
#' @importFrom methods setMethod
setMethod("[",
          signature(x = "cfDataList", i = "ANY", j = "ANY", drop = "ANY"),
          function (x, i, j){
            if (!missing(j))
              warning("column subscripts ignored")
            x@.Data[[i]]
          }
)

#' @importFrom methods setMethod
#' @rdname Extract
#' @aliases [[,cfDataList-method
setMethod("[[",
          signature(x = "cfDataList"),
          function (x, i){
            x@.Data[[i]]
          }
)
