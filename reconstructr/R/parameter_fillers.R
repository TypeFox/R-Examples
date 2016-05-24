#'@title automatically generate plausible padding values
#'@description \code{\link{padding_value}} is designed to automatically generate
#'acceptable values for the \code{padding_value} argument in \code{\link{session_length}}.
#'It does this by applying a statistical function to the range of inter-event values within
#'the dataset.
#'
#'@param sessions a list of sessions, generated with \code{\link{reconstruct_sessions}}
#'
#'@param method the method to use to generate the value. Options are "geometric mean",
#'"arithmetic mean", "median" or "other" (see below).
#'
#'@param fun the function to apply to generate \code{\link{padding_value}}, in the
#'event that you choose "other" for \code{method}.
#'
#'@param ... additional arguments to pass to fun, or the function called by the
#'method.
#'
#'@examples
#'
#'#Generate a padding value based on the geometric mean of the known times between events
#'data("session_dataset")
#'session_dataset$timestamp <- to_seconds(x = session_dataset$timestamp, format = "%Y%m%d%H%M%S")
#'events_by_user <- split(session_dataset$timestamp, session_dataset$UUID)
#'sessions <- reconstruct_sessions(events_by_user)
#'padding_val <- padding_value(sessions, "geometric mean")
#'padding_val
#'#41.07547
#'
#'@export
padding_value <- function(sessions, method, fun, ...){
  
  #Grab intertimes
  intertimes <- unlist(c_time_per_event(sessions))
  intertimes <- intertimes[intertimes > -1]
  
  #Apply!
  fun_apply <- switch(method,
                      "arithmetic mean" = mean,
                      "geometric mean" = function(x){exp(sum(log(x[x > 0])) / length(x))},
                      "median" = median,
                      fun)
  
  result <- fun_apply(intertimes, ...)
  return(result)
}