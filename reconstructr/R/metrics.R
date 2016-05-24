#'@title calculate the bounce rate within a session dataset
#'@description calculates the "bounce rate" within a set of sessions - the proportion of sessions
#'consisting only of a single event.
#'
#'@param sessions a list of sessions, generated with \code{\link{reconstruct_sessions}}
#'
#'@param decimal_places the number of decimal places to round the output to - set to 2 by default.
#'
#'@return a single numeric value, representing the percentage of sessions that are bounces.
#'
#'@seealso \code{\link{session_events}} for generaliseable event-level calculations, and
#'\code{link{event_time}} for performing operations on the time between events.
#'
#'@examples
#'#Calculate the bounce rate in the provided dataset.
#'#Load, convert timestamps to seconds, split
#'data("session_dataset")
#'session_dataset$timestamp <- to_seconds(x = session_dataset$timestamp, format = "%Y%m%d%H%M%S")
#'events_by_user <- split(session_dataset$timestamp, session_dataset$UUID)
#'
#'#Sessionise and calculate bounce rate
#'sessions <- reconstruct_sessions(events_by_user)
#'bounce_rate(sessions)
#'#[1]58
#'@export
bounce_rate <- function(sessions, decimal_places = 2){
  events <- session_events(sessions)
  return(round((length(events[events == 1])/(length(events))), digits = decimal_places)*100)
}

#'@title calculate the time between each event in a session, or set of sessions
#'@description As well as \code{\link{session_length}}, which calculates the length
#'of a session, \code{event_time} extracts the length of time between each event
#'(which can be used for purposes such as counting the likely time-on-page, in the
#'context of web analytics). It allows you to customise
#'the output format and optionally run analytical functions against
#'the results before they are returned
#'
#'@param sessions a list of sessions, generated with \code{\link{reconstruct_sessions}}
#'
#'@param as_vector whether to unlist the results and return them as a vector,
#'or return them as a list. Set to TRUE (vector) by default. In the event that you
#'are setting a value for "fun", your choice may impact the way the function
#'operates over the results
#'
#'@param fun an optional parameter indicating a function to pass over the results,
#'before they are returned.
#'
#'@param ... optional arguments to pass to fun
#'
#'@return a list, a vector, or the output format of \code{fun}
#'
#'@seealso \code{\link{session_length}} for calculating the length of time
#'spent within a session
#'
#'@examples
#'#Load, convert timestamps to seconds, split
#'data("session_dataset")
#'session_dataset$timestamp <- to_seconds(x = session_dataset$timestamp, format = "%Y%m%d%H%M%S")
#'events_by_user <- split(session_dataset$timestamp, session_dataset$UUID)
#'sessions <- reconstruct_sessions(events_by_user)
#'
#'#Extract the time between events, separating out the results for each session
#'list_event_time <- event_time(sessions, as_vector = FALSE)
#'
#'#Extract the arithmetic mean time between events
#'mean_event_time <- event_time(sessions, as.vector = TRUE, fun = mean, trim = 0.5)
#'
#'@export
event_time <- function(sessions, as_vector = TRUE, fun, ...){
  
  #Grab intertimes
  intertimes <- c_time_per_event(sessions)
  
  if(as_vector){
    intertimes <- unlist(intertimes)
  }
  
  #And now we play the control flow game!
  if(missing("fun")){
    return(intertimes)
  } else {
    return(fun(sessions, ...))
  }
}