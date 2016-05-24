#' switch.transition Add a transition parameter on a set of parameters or remove it
#' @title Add a transition parameter on a set of parameters or remove it
#' @author Marc Girondot
#' @return A vector with parameters
#' @param parameters A vector with parameters
#' @description Add a transition parameter on a set of parameters or remove it
#' @examples
#' \dontrun{
#' data(resultNest_6p)
#' # Get a set of parameters without transition
#' x1 <- resultNest_6p$par
#' # Generate a set of parameters with transition
#' x2 <- switch.transition(x1)
#' # Generate a set of parameters without transition
#' x3 <- switch.transition(x3)
#' }
#' @export

switch.transition <- function(parameters=stop("A set of parameters must be supplied")) {

if (any(substr(names(parameters), nchar(names(parameters))-1, nchar(names(parameters)))=="_L")) {
	xp <- parameters[substr(names(parameters), nchar(names(parameters))-1, nchar(names(parameters)))!="_L"]
	xp <- xp[substr(names(xp), 1, 10)!="transition"]
} else {
	xp <- parameters
	names(xp) <- paste0(names(xp), "_L")
	xp <- c(parameters, xp, transition_P=20, transition_S=0.01)
}


return(xp)
}

