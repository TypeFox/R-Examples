#' Event Range
#' 
#' Returns the minimum and maximum endpoints of all the events in an event table.
#' 
#' @param e An event table.
#' @seealso \code{\link{event_coverage}} for an alternative that accounts for gaps.
#' @export
#' @examples
#' event_range(events(1:5))            # no gaps
#' event_range(events(c(1,5), c(1,5))) # gaps
event_range <- function(e) {
  return(as_events(range(e$from, e$to)))
}

#' Event Coverage
#' 
#' Returns the intervals over which the number of events is always one or greater.
#' 
#'  If \code{closed = TRUE}, breaks between adjacent events are dropped. If \code{closed = FALSE}, breaks between adjacent events are retained, including point events on line event endpoints. Duplicate points are dropped in both cases.
#'
#' @param e An event table.
#' @param closed Logical value indicating whether events should be interpreted as closed intervals. If \code{TRUE}, coverage is continuous at breaks between two adjacent events.
#' @seealso \code{\link{event_gaps}} for gaps (the inverse of coverage), \code{\link{event_range}} for range (coverage with gaps ignored).
#' @export
#' @examples
#' e <- events(c(1, 2, 4, 8), c(3, 4, 5, 10))
#' event_coverage(e, closed = TRUE)  # retains breaks
#' event_coverage(e, closed = FALSE) # drops breaks
#' e <- events(c(0, 2, 2, 2, 8, 10), c(0, 2, 2, 6, 10, 10))
#' event_coverage(e, closed = TRUE)  # retains isolated points
#' event_coverage(e, closed = FALSE) # retains isolated points and points adjacent to lines
event_coverage <- function(e, closed = TRUE) {
  e.range <- event_range(e)
  e.gaps <- event_gaps(e, closed = closed)
  e.coverage <- event_gaps(e.gaps, closed = FALSE, range = e.range)
  # restore isolated boundary point events
  add <- e.range != event_range(e.coverage)
  # restore adjacent boundary point events
  if (!closed) {
    add[1] <- add[1] || any(rowSums(e[c("from", "to")] == e.range[[1]]) == 2)
    add[2] <- add[2] || any(rowSums(e[c("from", "to")] == e.range[[2]]) == 2)
  }
  if (any(add)) {
    e.coverage <- rbind(if_else(add[[1]], rep(e.range[[1]], 2), NULL), e.coverage, if_else(add[[2]], rep(e.range[[2]], 2), NULL))
  }
  return(e.coverage)
}

#' Event Gaps
#' 
#' Returns the intervals over which there are no events.
#'
#' @param e An event table.
#' @param closed Logical value indicating whether events should be interpreted as closed intervals. If \code{TRUE}, no gaps are returned at breaks between two adjacent events.
#' @param range An event table specifying, by its \code{\link{event_range}}, the interval within which to check for gaps. If \code{NULL}, the range of \code{e} is used.
#' @seealso \code{\link{event_coverage}} for coverage (the inverse of gaps), \code{\link{fill_event_gaps}} for filling gaps with empty events.
#' @export
#' @examples
#' event_gaps(events(c(1, 3, 5), c(2, 4, 5)))    # gaps between events
#' event_gaps(events(1:5))                       # no gaps
#' event_gaps(events(1:5), closed = FALSE)       # gaps at breaks
#' event_gaps(events(1:5), range = events(0, 6)) # gaps to edge of range  
event_gaps <- function(e, closed = TRUE, range = NULL) {
  # Crop to bounds if not event range
  if (!is.null(range)) {
    range <- event_range(range)
    if (!nrow(e)) {
      # no data in range
      return(range)
    }
    e.range <- event_range(e)
    if (!(range$from <= e.range$from && range$to >= e.range$to)) {
      # bounds intersect event range
      e <- crop_events(e, range)
    }
    if (range$from < e.range$from) {
      # bounds extend past event range (from)
      e <- rbind(rep(range$from, 2), e[c("from", "to")])
    }
    if (range$to > e.range$to) {
      # bounds extend past event range (to)
      e <- rbind(e[c("from", "to")], rep(range$to, 2))
    }
  }
  # Track overlaps by extending events by cumulative max
  if (is_unsorted_events(e)) {
    e <- sort_events(e)
  }
  e$to <- cummax(e$to)
  # Gaps occur when from[i+1] > to[i] if closed intervals
  # (and when == if open intervals)
  if (closed) { 
    isgap <- which(e$from[-1] > e$to[-nrow(e)])
  } else { 
    isgap <- which(e$from[-1] >= e$to[-nrow(e)])
  }
  # Build gaps event table
  # (remove possible duplicate point gaps)
  return(unique(events(e$to[isgap], e$from[isgap + 1])))
}

#' Event Overlaps
#' 
#' Returns the number of events on each interval. Useful for sampling the original data with \code{\link{sample_events}} at the highest possible resolution that nevertheless flattens overlapping events.
#' 
#' Point events are preserved and line events are cut as necessary at the endpoints of other point or line events.
#' 
#' @param e An event table.
#' @return An endpoint-only event table with column "n" listing the number of overlapping events on that interval.
#' @seealso \code{\link{event_coverage}}.
#' @export
#' @examples
#' e <- events(c(0, 10, 15, 25, 30), c(10, 20, 25, 40, 30))
#' event_overlaps(e)
event_overlaps <- function(e) {
  # Sort event table
  if (is_unsorted_events(e)) {
    e <- sort_events(e)
  }
  # Cut events at event endpoints
  # (don't cut new points out points, just cut lines)
  e.cut <- cut_events(e, c(e$from, e$to))
  temp <- stats::aggregate(e.cut$from, by = e.cut[c("from","to")], FUN = length)
  names(temp)[3] <- "n"
  return(temp)
}