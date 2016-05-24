#' Subset results to those within specified area around zip code.
#'
#' @param sccall Current list of parameters carried forward from prior
#'     functions in the chain (ignore)
#' @param zip A zipcode
#' @param distance An integer distance in miles or kilometers
#' @param km A boolean value set to \code{TRUE} if distance should be
#'     in kilometers (default is \code{FALSE} for miles)
#' @examples
#' \dontrun{
#' sc_zip(37203)
#' sc_zip(37203, 50)
#' sc_zip(37203, 50, km = TRUE)
#' }

#' @export
sc_zip <- function(sccall, zip, distance = 25, km = FALSE) {

    ## check first argument
    if (identical(class(try(sccall, silent = TRUE)), 'try-error')
        || !is.list(sccall)) {
        stop('Chain not properly initialized. Be sure to start with sc_init().',
             call. = FALSE)
    }

    ## check second argument
    if (missing(zip) || !is.numeric(zip) || nchar(zip) != 5) {
        stop('Must provide a 5-digit zip code.', call. = FALSE)
    }

    stub <- '&_zip=' %+% zip %+% '&_distance=' %+% distance

    if (km) {
        sccall[['zip']] <- stub %+% 'km'
    } else {
        sccall[['zip']] <- stub
    }

    sccall

}



