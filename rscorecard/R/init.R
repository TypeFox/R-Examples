#' Initialize chained request.
#'
#' This function initializes the data request. It should always be
#' the first in the series of piped functions.
#'
#' @param dfvars Set to \code{TRUE} if you would rather use the
#'     developer-friendly variable names used in actual API call.
#'
#' @examples
#' \dontrun{
#' sc_init()
#' sc_init(dfvars = TRUE)
#' }

#' @export
sc_init <- function(dfvars = FALSE) {

    list('dfvars' = dfvars,
         'select' = NULL,
         'filter' = NULL,
         'zip' = NULL,
         'year' = 2013)

}
