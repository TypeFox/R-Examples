#############################################################################
# style.R
#############################################################################

#' @include style_apa.R
NULL

#########
# initialise function
#########

#' Setting the publication style functions
#'
#' A list is returned which contains all the default functions for the
#' specified publication style.
#'
#' @param x a character indicating the publication style.
#'
#' @seealso See \code{\link{pp_opts_style}} to get and set all currently used
#' style functions.
#'
#' @examples
#' pp_opts_style$set(pp_init_style())
#' pp_opts_style$set(pp_init_style("apa"))
#' 
#' @export
pp_init_style <- function(x = c("apa"))
{
    x <- match.arg(x)

    ret <- switch(EXPR = x,
                  "apa" = style.apa.init())

    return(ret)
}

#########
# calling functions
#########
