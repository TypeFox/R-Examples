#' Fill Event Gaps
#' 
#' \code{fill_event_gaps} fills gaps below a maximum length with empty events. \code{collapse_event_gaps} shifts event endpoints to close gaps below a maximum length.
#' 
#' @param e An event table.
#' @param max.length The maximum length of gaps to be filled or closed.
#' @seealso \code{\link{event_gaps}}
#' @export
#' @rdname fill_event_gaps
#' @examples
#' e <- events(c(1, 4), c(2, 5), x = 1)
#' fill_event_gaps(e)
#' fill_event_gaps(e, max.length = 1)
#' collapse_event_gaps(e)
#' collapse_event_gaps(e, max.length = 1)
fill_event_gaps <- function(e, max.length = Inf) {
  e.gaps <- event_gaps(e)
  gaps.length <- e.gaps$to - e.gaps$from
  fill <- gaps.length < max.length
  ngaps <- sum(fill)
  nevents <- nrow(e)
  if (ngaps > 0) {
    e <- e[c(seq_len(nrow(e)), rep(NA, ngaps)), , drop = FALSE]
    e[nevents + (1:ngaps), c("from", "to")] <- e.gaps[fill, ]
  }
  return(e)
}
#' @export 
#' @rdname fill_event_gaps
collapse_event_gaps <- function(e, max.length = Inf) {
  if (is_unsorted_events(e)) {
    ids <- order(e$from, e$to)
    e <- e[ids, , drop = FALSE]
    reorder <- TRUE
  } else {
    reorder <- FALSE
  }
  ne <- nrow(e)
  cummax.to <- cummax(e$to)
  gap <- e$from[-1] - cummax.to[-ne]
  is.gap <- which(gap > 0 & gap < max.length)
  temp <- numeric(ne)
  temp[is.gap + 1] <- gap[is.gap]
  e[c("from", "to")] <- e[c("from", "to")] - cumsum(temp)
  # Reorder as needed
  if (reorder) {
    return(e[order(ids), , drop = FALSE])
  } else {
    return(e)
  }
}

#' Transform Events
#' 
#' Transforms events by scaling, then translating their endpoint positions. That is, the transformed \code{[from, to] = scale * [from, to] + translate}.
#' 
#' @param e An event table.
#' @param scale Number by which event endpoints should be scaled.
#' @param translate Number by which event endpoints should be translated.
#' @export
#' @examples
#' e <- events(c(10, 100), c(100, 1000))
#' transform_events(e, scale = 2, translate = 1)
transform_events <- function(e, scale = 1, translate = 0) {
  e[c("from","to")] <- scale[1] * e[c("from","to")] + translate[1]
  return(e)
}

#' Crop Events
#' 
#' Crops events to the specified intervals. Events are cut at interval endpoints and any whole or partial events lying outside the intervals are removed.
#'
#' @param e An event table.
#' @param crops An event table specifying the intervals for cropping. Point intervals are allowed, and will create new point events where they intersect the interior, but not the endpoints, of line events.
#' @param scaled.cols Names or indices of the columns of the event table to be rescaled after cutting (see \code{\link{cut_events}}). Names are interpreted as regular expressions (\code{\link{regex}}) matching full column names.
#' @seealso \code{\link{cut_events}} for only cutting events.
#' @export
#' @examples
#' e <- events(c(0, 10, 20), c(10, 20, 30), x = 10)
#' crop_events(e, events(c(0, 15)))
#' crop_events(e, events(c(0, 5, 15)))
#' crop_events(e, events(c(0, 5, 15)), scaled.cols = "x")
#' crop_events(e, events(c(0, 5, 5, 15)), scaled.cols = "x")   # creates new points inside lines
#' crop_events(e, events(c(0, 10, 10, 15)), scaled.cols = "x") # but not at line event endpoints
crop_events <- function(e, crops, scaled.cols = NULL) {  
  e.cut = cut_events(e, cuts = crops, scaled.cols = scaled.cols)
  inx = find_intersecting_events(crops, e.cut)
  keep = .rowSums(inx, nrow(inx), ncol(inx)) > 0
  return(e.cut[keep, , drop = FALSE])
}

#' Cut Events
#' 
#' Cuts events at the specified locations.
#'
#' Line events straddling cut locations are cut into multiple event segments. Columns \code{scaled.cols} are scaled by the fraction of the original event length in each resulting event (which assumes that these variables were uniformly distributed over the original interval). To have a record of the parents of the resulting event segments, append an unique identification field to the event table before calling this function.
#'
#' @param e an event table.
#' @param cuts the cut locations. Can be either a numeric vector or an event table. If an event table that contains points, point events will be created where they intersect the interior, but not the endpoints, of line events in \code{e}.
#' @param scaled.cols names or indices of the event table columns to be scaled to their new length after cutting. Names are interpreted as regular expressions (\code{\link{regex}}) matching full column names.
#' @seealso \code{\link{crop_events}} for both cutting and removing events.
#' @export
#' @examples
#' e <- events(c(0, 10, 20), c(10, 20, 30), x = 10)
#' cut_events(e, events(c(0, 5, 15)))
#' cut_events(e, events(c(0, 5, 15)), scaled.cols = "x")
#' cut_events(e, events(c(0, 5, 5, 15)), scaled.cols = "x")   # creates new points inside lines
#' cut_events(e, events(c(0, 10, 10, 15)), scaled.cols = "x") # but not at line event endpoints
cut_events <- function(e, cuts, scaled.cols = NULL) {
  # Initialize inputs
  if (is.numeric(cuts)) {
    cuts <- sort(unique(cuts))
  } else {
    # remove duplicate intervals
    cuts <- unique(cuts[c("from", "to")])
    # for lines, unlist and remove duplicates
    # for points, unlist and combine
    is.lines <- cuts$to != cuts$from
    cuts <- sort(c(unique(unlist(cuts[is.lines, ])), unlist(cuts[!is.lines, ])))
    # limit repeating cutpoint sequence to 2 (a point).
    reps <- rle(cuts)
    reps$lengths[reps$lengths > 2] = 2
    cuts <- rep(reps$value, reps$lengths)
  }
  if (is.character(scaled.cols)) {
    scaled.cols <- unique(unlist(rgrep_exact(scaled.cols, names(e))))
  }
  if (!length(scaled.cols)) {
    scaled.cols <- NULL
  }
  # Find events straddling the cuts
  # (leave point events out of this)
  incuts <- find_intersecting_events(events(cuts, cuts), e, equal.points = FALSE)
  ncuts <- rowSums(incuts)
  if (sum(ncuts)) {
    # Expand events to accomodate new event segments
    ei <- seq_len(nrow(e))
    e <- e[rep(ei, ncuts + 1), ]
    # Indices of cut segments
    ind <- ncuts > 0
    ncuts <- ncuts[ind]
    i <- ei[ind]
    i0 <- i + c(0, cumsum(ncuts))[seq_along(i)]
    i1 <- i0 + ncuts
    # For each cut event
    for (n in seq_along(i)) {
      # Update measures
      pts <- c(e$from[i0[n]], cuts[incuts[i[n],]], e$to[i0[n]])
      e$from[i0[n]:i1[n]] <- pts[-length(pts)]
      e$to[i0[n]:i1[n]] <- pts[-1]
      # Rescale columns
      if (!is.null(scaled.cols)) {
        l <- diff(pts)
        e[i0[n]:i1[n], scaled.cols] = e[i0[n]:i1[n], scaled.cols] * (l / sum(l))
      }
    }
  }
  return(e)
}