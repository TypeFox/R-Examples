#' Select scorecard data year.
#'
#' This function is used to select the year of the data.
#'
#' @param sccall Current list of parameters carried forward from prior
#'     functions in the chain (ignore)
#' @param year Four-digit year (default is 2013)
#'
#' @section Important notes:
#' \enumerate{
#' \item Not all variables have a year option.
#' \item At this time, only one year at a time is allowed.
#' \item The year selected is not necessarily the year the data were produced.
#' It may be the year the data were collected. For data collected over split
#' years (fall to spring), it is likely the year represents the fall data (\emph{e.g.,}
#' 2011 for 2011/2012 data).
#' }
#'
#' Be sure to check with the College Scorecard
#' \href{https://collegescorecard.ed.gov/assets/FullDataDocumentation.pdf}{data
#' documentation report} when choosing the year.
#'
#' @examples
#' \dontrun{
#' sc_year()
#' sc_year(2012)
#' }

#' @export
sc_year <- function(sccall, year) {

    ## check first argument
    if (identical(class(try(sccall, silent = TRUE)), 'try-error')
        || !is.list(sccall)) {
        stop('Chain not properly initialized. Be sure to start with sc_init().',
             call. = FALSE)
    }

    ## check second argument
    if (missing(year) || !is.numeric(year) || year < 1900 || year > 2099) {
        stop('Must provide a 4-digit year in 1900s or 2000s.', call. = FALSE)
    }

    ## get vars
    sccall[['year']] <- year
    sccall

}


