#' sample 'msl.trend' object
#'
#' Output of call to \code{\link{msl.trend}} used extensively in examples throughout
#' this Manual.
#'
#' @docType data
#'
#' @usage data(s)
#'
#' @format msl.trend object
#'
#' @details This \code{\link{msl.trend}} object is used extensively in the
#' examples throughout this manual in order to call the object direct rather than
#' producing the same via original code which can be computationally expensive. This
#' object results from a decomposition of the Baltimore record, filling gaps with
#' spline interpolation and using 500 iterations to generate error margins via
#' bootstrapping.
#'
#' \strong{Note: }Ordinarily the user would call 'File.csv' direct from working
#' directory, creating the 'msl.trend' object using the following sample code:
#'
#' s <- msl.trend('Balt.csv', fillgaps = 3, iter = 500, 'BALTIMORE, USA') # DON'T RUN
#'
#' @seealso \code{\link{msl.trend}}, \code{\link{msl.forecast}},
#' \code{\link{msl.plot}}, \code{\link{msl.pdf}}, \code{\link{summary}},
#' \code{\link{Balt}}.
#'
#' @examples
#' data(s)
#' str(s) # check structure of object
"s"
