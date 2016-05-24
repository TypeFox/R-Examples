#' initiates databel object
#'
#' this is a simple wrapper for the "new" function
#' creating a databel object
#'
#' @param baseobject name of the file or \link{databel-class} object
#' @param cachesizeMb cache size (amount of RAM) to be used
#' @param readonly readonly flag
#'
#' @author Yurii Aulchenko
#' @export
#'

databel <- function(baseobject, cachesizeMb=64, readonly=TRUE)
{
#	if (missing(cachesizeMb)) {
#		if (is(baseobject,"databel"))
#		{
#			cachesizeMb <- cachesizeMb(baseobject)
#		} else if (is(baseobject,"character")) {
#			cachesizeMb <- 64
#		}
#	}
    ret <- new(Class="databel", baseobject=baseobject,
               cachesizeMb=cachesizeMb, readonly=readonly);
    return(ret)
}
