#' Split Data Frame into Multiple Groups
#'
#' A data frame is split into a list of data subsets defined by multiple groups.
#'
#' @param data data frame, matrix or vector containing several subsets of data
#' @param grouping either a character vector containing the names of the grouping variables
#'   or a model formula specifying dependent,
#'   independent and grouping variables in the form:
#'   \code{dependent ~ independent | group1 + group2 + ...}.
#'   It may also be a factor or list of factors as in \code{\link{split}}.
#' @param drop if drop is TRUE, unused factor levels are dropped from the result.
#'   The default is to drop all factor levels.
#' @param sep	string to construct the new level labels by joining the
#'   constituent ones.
#' @param \dots other parameters passed to \code{\link{split}}, see details.
#'
#' @details This function is wrapper around \code{\link{split}} with
#'   different defaults, slightly different behavior, and methods for additional
#'   argument classes. \code{multisplit} returns always a data frame.
#'
#' @return list containing data frames of the data subsets as its elements.
#'   The components of the list are named by their grouping levels.
#'
#' @seealso \code{\link{split}}
#'
#' @examples
#'
#'
#'data(bactgrowth)
#'
#'## simple method
#'spl <- multisplit(bactgrowth, c("strain", "conc", "replicate"))
#'
#'## preferred method
#'spl <- multisplit(bactgrowth, value ~ time | strain + conc + replicate)
#'
#'## show what is in one data set
#'spl[[1]]
#'summary(spl[[1]])
#'
#'## use factor combination
#'spl[["D:0:1"]]
#'summary(spl[["D:0:1"]])
#'
#'
#'lapply(spl, FUN=function(x)
#'  plot(x$time, x$value,
#'       main=paste(x[1, "strain"], x[1, "conc"], x[1, "replicate"], sep=":")))
#'
#'
#' @rdname multisplit
#' @keywords internal
#' @exportMethod multisplit
#'
setMethod("multisplit", c("data.frame", "formula"),
          function(data, grouping, drop = TRUE, sep=":", ...) {
            if (missing(grouping) || (length(grouping) != 3L))
              stop("'grouping' missing or incorrect")

            if (is.matrix(data))
              data <- as.data.frame(data)
            if (!is.data.frame(data))
              stop("'data' must be a data frame or matrix")

            p <- parse_formula(grouping)
            nm <- names(data)
            if (is.na(p$groups[1]))
              stop("grouping variable(s) missing")
            if (!all(c(p$valuevar, p$timevar) %in% names(data)))
              stop("dependent and independent variables must be column names of data")
            if (!all(p$groups %in% names(data)))
              stop("all grouping criteria must be column names of data")

            split(data[c(p$timevar, p$valuevar)], data[p$groups],
                  drop = drop, sep = sep)
          })


#' @rdname multisplit
#' @keywords internal
#' @exportMethod multisplit
#'
setMethod("multisplit", c("data.frame", "character"),
          function(data, grouping, drop = TRUE, sep = ":", ...) {
            if (!all(grouping %in% names(data)))
              stop("Not all groups found in data frame")
            split(data, data[, grouping], drop = drop, sep = sep)
          })

#' @rdname multisplit
#' @keywords internal
#' @exportMethod multisplit
#'
setMethod("multisplit", c("data.frame", "factor"),
          function(data, grouping, drop = TRUE, sep = ":", ...) {
            split(data, list(grouping), drop = drop, sep = sep, ...)
          })

#' @rdname multisplit
#' @keywords internal
#' @exportMethod multisplit
#'
setMethod("multisplit", c("data.frame", "list"),
          function(data, grouping, drop = TRUE, sep = ":", ...) {
            split(data, grouping, drop = drop, sep = sep, ...)
          })

#' @rdname multisplit
#' @keywords internal
#' @exportMethod multisplit
#'
setMethod("multisplit", c("ANY", "ANY"),
          function(data, grouping, drop = TRUE, sep = ":", ...) {
              data <- as.data.frame(data)
              multisplit(data, grouping, drop = drop, sep = sep, ...)
            })

