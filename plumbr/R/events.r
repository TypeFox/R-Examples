#' Plumbr events
#'
#' Plumbr data structures send only single event for data changes:
#' data_changed.  This has a two arguments, i and j.  Either both are NULL,
#' indicating a change in the shape of the underlying data, or they give the
#' the locations of changed data values.
#'
#' @param mf muta frame
#' @param callback function with arguments i and j
#' @export
add_listener <- function(mf, callback) {
  changed(mf)$connect(callback)
}

##' Remove a listener, identified by the ID returned by
##' \code{\link{add_listener}}. 
##'
##' @param mf mutaframe
##' @param id value returned by \code{add_listener} when originally
##' connecting the handler
##' @export
remove_listener <- function(mf, id) {
  changed(mf)$disconnect(id)
}

#' Is the event a shape changed event?
#' 
#' @param i col index
#' @param j row index
#' @export
shape_changed <- function(i, j) is.null(i) || is.null(j)

#' Combine list of events into single event.
#'
#' If any event is a \code{shape_changed} event, return it.  Otherwise, 
#' take the unique elements of the union of all element changes.
#'
#' @param events a list of event parameters
#' @return a unified event
combine_data_events <- function(events) {
  for(event in events) {
    if (shape_changed(event$i, event$j)) return(list(i = NULL, j = NULL))
  }
  
  unique(do.call("rbind", lapply(events, as.data.frame, 
    stringsAsFactors = FALSE)))
}

#' Get the 'changed' signal
#' @param mf a mutaframe
changed <- function(mf) {
  attr(mf, "changed")
}

#' Notify listeners that data has changed.
#' @param mf mutaframe
#' @param i,j row and column indices
#' @export
notify_listeners <- function(mf, i, j) {
  changed(mf)$emit(i, j)
}

#' Pause (cache) events.
#' 
#' When a mutaframe is paused, it accumulates events without passing them on.
#' When unpaused, it accumulates all events into a single event and passes it
#' on.
#'
#' This is a performance optimisation for when you expect many changes: 
#' pause the mutaframe, perform all the changes and then unpause.
#'
#' @param mf mutaframe
#' @export
pause <- function(mf) {
  changed(mf)$buffer()
}

#' Unpause (reply) events.
#'
#' @param mf mutaframe
#' @export
unpause <- function(mf) {
  changed(mf)$flush()
}

#' Is a mutaframe currently paused?
#' @param mf a mutaframe
#' @export
is_paused <- function(mf) changed(mf)$paused
