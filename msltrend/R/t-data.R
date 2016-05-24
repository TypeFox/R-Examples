#' sample 'msl.forecast' object
#'
#' Output of call to \code{\link{msl.forecast}} used extensively in examples throughout
#' this Manual.
#'
#' @docType data
#'
#' @usage data(t)
#'
#' @format msl.forecast object
#'
#' @details This \code{\link{msl.forecast}} object is used extensively in the
#' examples throughout this manual in order to call the object direct rather than
#' producing the same via original code which can be computationally expensive. This
#' object results from a decomposition of the Baltimore record, filling gaps with
#' spline interpolation and using 500 iterations to generate error margins via
#' bootstrapping (see \code{\link{s}}). This 'msl.trend' object is then parsed to
#' \code{\link{msl.forecast}} with the addition of 1000 millimetres of sea level rise
#' between the end of the historical record and 2100.
#'
#' \strong{Note: }Ordinarily the user would call 'File.csv' direct from working
#' directory, creating the 'msl.trend' object first, then creating the above-mentioned
#' \code{\link{msl.forecast}} object using the following sample code:
#'
#' s <- msl.trend('Balt.csv', fillgaps = 3, iter = 500, 'BALTIMORE, USA') # DON'T RUN
#'
#' t <- msl.forecast(s, slr = 1000) # DON'T RUN
#'
#' @seealso \code{\link{msl.trend}}, \code{\link{msl.forecast}},
#' \code{\link{msl.plot}}, \code{\link{msl.pdf}}, \code{\link{summary}},
#' \code{\link{Balt}}, \code{\link{s}}.
#'
#' @examples
#' data(t)
#' str(t) # check structure of object
"t"
