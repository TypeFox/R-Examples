
######################################################################################################################

# Function: Event.
# Argument: A list or vector of numeric.
# Description: This function is used to create an object of class Event.
#' @export
Event = function(n.events, rando.ratio=NULL) {

  # Error checks
  if (any(!is.numeric(unlist(n.events)))) stop("Event: number of events must be numeric.")
  if (any(unlist(n.events) %% 1 !=0)) stop("Event: number of events must be integer.")
  if (any(unlist(n.events) <=0)) stop("Event: number of events must be strictly positive.")

  if (!is.null(rando.ratio)){
    if (any(!is.numeric(unlist(rando.ratio)))) stop("Event: randomization ratio must be numeric.")
    if (any(unlist(rando.ratio) %% 1 !=0)) stop("Event: randomization ratio must be integer.")
    if (any(unlist(rando.ratio) <=0)) stop("Event: randomization ratio must be strictly positive.")
  }

  event = list(n.events = unlist(n.events), rando.ratio = unlist(rando.ratio))

  class(event) = "Event"
  return(event)
  invisible(event)

}